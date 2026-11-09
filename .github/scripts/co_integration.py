"""Run a Code Ocean integration pipeline at a commit and report its verdict.

Contract with the pipeline:

- It accepts the commit under test as a named parameter (``--commit-param``).
- Its final step writes ``integration_verdict.json`` at the root of its results::

      {
        "schema_version": 1,
        "status": "pass" | "fail" | "error",
        "commit": "<commit that ran>",
        "summary": "<optional one line for the pull request>",
        "datasets": [{"name": "...", "status": "pass" | "fail" | "error", "message": "..."}]
      }

  ``fail`` means it ran and a check failed; ``error`` means something could not run.
- It exits nonzero when it cannot finish. It may exit 0 with a ``fail`` verdict, since a
  failing final step can stop Nextflow from publishing the verdict.

The check passes when the pipeline exits 0, the verdict status is ``pass``, and the
verdict commit matches the commit under test. A missing or unreadable verdict is an error.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

from codeocean import CodeOcean
from codeocean.computation import ComputationState, NamedRunParam, RunParams
from codeocean.error import Error as CodeOceanError

CO_DOMAIN = "https://codeocean.allenneuraldynamics.org"
POLL_SECONDS = 60
LOG_FILE = "output"
VERDICT_FILE = "integration_verdict.json"
SCHEMA_VERSION = 1
PASS = "pass"
LOG_TAIL_LINES = 40


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pipeline-id", required=True, help="Code Ocean pipeline id.")
    parser.add_argument("--commit", required=True, help="Commit under test.")
    parser.add_argument("--commit-param", default="test_commit", help="Pipeline parameter that receives the commit.")
    parser.add_argument("--timeout-minutes", type=float, default=330, help="Maximum run time to wait.")
    return parser.parse_args()


def read_result(client: CodeOcean, computation_id: str, path: str) -> str | None:
    """Return the text of a result file, or None when it cannot be read."""
    try:
        url = client.computations.get_result_file_urls(computation_id, path).download_url
        with urllib.request.urlopen(url) as response:
            return response.read().decode()
    except (CodeOceanError, urllib.error.URLError):
        return None


def start_run(client: CodeOcean, args: argparse.Namespace):
    """Start the pipeline and confirm Code Ocean recorded the commit parameter."""
    computation = client.computations.run_capsule(
        RunParams(
            pipeline_id=args.pipeline_id,
            named_parameters=[NamedRunParam(param_name=args.commit_param, value=args.commit)],
        )
    )
    recorded = {p.param_name: p.value for p in computation.parameters or [] if p.param_name}
    if recorded.get(args.commit_param) != args.commit:
        raise SystemExit(
            f"Code Ocean recorded {args.commit_param}={recorded.get(args.commit_param)!r}, expected "
            f"{args.commit!r}; check the parameter name in the pipeline's app panel."
        )
    return computation


def run_failures(computation) -> list[str]:
    """Return why the run itself failed; empty when it completed with exit code 0."""
    if computation.state != ComputationState.Completed:
        return [f"computation state is {computation.state.value}"]
    if computation.exit_code != 0:
        return [f"pipeline exited with code {computation.exit_code}"]
    return []


def parse_verdict(text: str | None) -> tuple[dict | None, list[str]]:
    """Return the verdict and any contract violations in it."""
    if text is None:
        return None, [f"{VERDICT_FILE} is missing"]
    try:
        verdict = json.loads(text)
    except json.JSONDecodeError as error:
        return None, [f"{VERDICT_FILE} is not valid JSON: {error}"]
    if not isinstance(verdict, dict):
        return None, [f"{VERDICT_FILE} is not a JSON object"]
    if verdict.get("schema_version") != SCHEMA_VERSION:
        return verdict, [f"{VERDICT_FILE} schema_version is {verdict.get('schema_version')!r}, expected {SCHEMA_VERSION}"]
    return verdict, []


def verdict_failures(verdict: dict, commit: str) -> list[str]:
    """Return failures stated by or implied by the verdict."""
    failures = []
    if verdict.get("status") != PASS:
        failures.append(f"verdict status is {verdict.get('status')!r}")
    if verdict.get("commit") != commit:
        failures.append(f"pipeline ran commit {verdict.get('commit')!r}, expected {commit!r}")
    return failures


def dataset_table(datasets: list) -> list[str]:
    """Return a markdown table of per-dataset results."""
    rows = [d for d in datasets if isinstance(d, dict)]
    if not rows:
        return []
    lines = ["", "| Dataset | Status | Message |", "| --- | --- | --- |"]
    for d in rows:
        message = str(d.get("message", "")).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {d.get('name', '?')} | {d.get('status', '?')} | {message} |")
    return lines


def write_step_summary(text: str) -> None:
    """Append text to the GitHub Actions job summary when running in Actions."""
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as f:
            f.write(text + "\n")


def main() -> int:
    """Run the pipeline, wait for it, report, and return a process exit code."""
    args = parse_args()
    client = CodeOcean(domain=CO_DOMAIN, token=os.environ["CODEOCEAN_TOKEN"])
    slug = client.capsules.get_capsule(args.pipeline_id).slug

    computation = start_run(client, args)
    run_url = f"{CO_DOMAIN}/capsule/{slug}/tree?results={computation.id}_computation"
    print(f"Started computation {computation.id} for {args.commit}: {run_url}", flush=True)

    verdict = None
    try:
        computation = client.computations.wait_until_completed(
            computation, polling_interval=POLL_SECONDS, timeout=args.timeout_minutes * 60
        )
    except TimeoutError:
        failures = [f"pipeline did not finish within {args.timeout_minutes} minutes"]
    else:
        failures = run_failures(computation)
        verdict, contract_failures = parse_verdict(read_result(client, computation.id, VERDICT_FILE))
        failures += contract_failures
        if verdict is not None and not contract_failures:
            failures += verdict_failures(verdict, args.commit)

    lines = [f"## Code Ocean integration check {'failed' if failures else 'passed'}"]
    if verdict and verdict.get("summary"):
        lines.append(str(verdict["summary"]))
    lines += [f"- Commit: `{args.commit}`", f"- Run: [{computation.id}]({run_url})"]
    lines += [f"- Failure: {f}" for f in failures]
    if verdict:
        lines += dataset_table(verdict.get("datasets", []))
    if failures:
        log = read_result(client, computation.id, LOG_FILE) or "log unavailable"
        lines += ["", "Pipeline log (tail):", "```", *log.splitlines()[-LOG_TAIL_LINES:], "```"]

    text = "\n".join(lines)
    print(text)
    write_step_summary(text)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

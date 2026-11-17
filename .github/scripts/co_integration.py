"""Run a Code Ocean integration capsule for each data asset at the given commit and report the outcomes.

Contract with the capsule:

- It accepts the commit under test as a named parameter (``--commit-param``).
- It exits 0 when every check passes and nonzero when a check fails or cannot run.

The check passes when every capsule run completes with exit code 0. A failed, cancelled, or timed-out run fails the check and cancels the runs still in flight.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

from aind_codeocean_dispatch import Job, Outcome, job_key, run_jobs
from codeocean import CodeOcean
from codeocean.computation import DataAssetsRunParam, NamedRunParam, RunParams
from codeocean.error import Error as CodeOceanError

CO_DOMAIN = "https://codeocean.allenneuraldynamics.org"
POLL_SECONDS = 60
LOG_FILE = "output"
LOG_TAIL_LINES = 40


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--capsule-id", required=True, help="Code Ocean Capsule id.")
    parser.add_argument("--commit", required=True, help="Commit under test.")
    parser.add_argument("--data-asset-ids", type=json.loads, default=[], help="JSON array of data asset ids.")
    parser.add_argument("--commit-param", default="test_commit", help="Capsule parameter that receives the commit.")
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


def start_runs(client: CodeOcean, args: argparse.Namespace) -> list[Outcome]:
    """Start one capsule run per data asset and wait for every outcome."""
    jobs = []
    for asset_id in args.data_asset_ids:
        params = RunParams(
            capsule_id=args.capsule_id,
            data_assets=[DataAssetsRunParam(id=asset_id, mount="data")],
            named_parameters=[NamedRunParam(param_name=args.commit_param, value=args.commit)],
        )
        jobs.append(Job(label=asset_id, run_params=params, key=job_key(params)))

    return run_jobs(
        client,
        jobs,
        launch_interval_s=1,
        start_retries=3,
        poll_interval_s=POLL_SECONDS,
        timeout_s=args.timeout_minutes * 60,
        cancel_on_failure=True,
        on_outcome=print,
    )


def write_step_summary(text: str) -> None:
    """Append text to the GitHub Actions job summary when running in Actions."""
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as f:
            f.write(text + "\n")


def main() -> int:
    """Run one capsule per data asset at the commit, report each outcome, and return a process exit code."""
    args = parse_args()
    client = CodeOcean(domain=CO_DOMAIN, token=os.environ["CODEOCEAN_TOKEN"])
    slug = client.capsules.get_capsule(args.capsule_id).slug

    outcomes = start_runs(client, args)

    tests_failed = False
    sections = []
    for outcome in outcomes:
        comp_id = outcome.computation_ids[-1]

        lines = [f"### {outcome.job.label}: {'passed' if outcome.ok else 'failed'}"]
        run_url = f"{CO_DOMAIN}/capsule/{slug}/tree?results={comp_id}_computation"
        lines.append(f"- Run: [{comp_id}]({run_url})")

        if not outcome.ok:
            tests_failed = True
            lines.append(f"- Failure: job {outcome.status}: {outcome.failure}: {outcome.error}")
            log = read_result(client, comp_id, LOG_FILE) or "log unavailable"
            lines += ["", "Capsule log (tail):", "```", *log.splitlines()[-LOG_TAIL_LINES:], "```"]

        sections.append("\n".join(lines))

    header = "\n".join(
        [f"## Code Ocean integration check {'failed' if tests_failed else 'passed'}", f"- Commit: `{args.commit}`"]
    )
    text = "\n\n".join([header, *sections])
    print(text)
    write_step_summary(text)
    return 1 if tests_failed else 0


if __name__ == "__main__":
    sys.exit(main())

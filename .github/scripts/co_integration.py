"""Run the GIAnT-MATLAB Code Ocean integration capsule at a commit and check its results."""

import argparse
import json
import os
import sys
import urllib.request

from codeocean import CodeOcean
from codeocean.computation import (
    ComputationEndStatus,
    ComputationState,
    DataAssetsRunParam,
    RunParams,
)

CO_DOMAIN = "https://codeocean.allenneuraldynamics.org"
POLL_SECONDS = 30
COMMIT_FILE = "giant_commit.txt"
TRIAL_TABLE_FILE = "trial_table.h5"
SUMMARY_FILE = "smoke_test.json"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--giant-sha", required=True, help="GIAnT-MATLAB commit to test.")
    parser.add_argument("--capsule-id", required=True, help="Integration capsule id.")
    parser.add_argument(
        "--session-asset-id",
        default=None,
        help="Raw SLAP2 session data asset id; omit to use the data attached to the capsule.",
    )
    parser.add_argument("--timeout-minutes", type=float, default=60, help="Maximum run time to wait.")
    return parser.parse_args()


def read_result_text(client: CodeOcean, computation_id: str, path: str) -> str:
    """Return the text of one result file."""
    url = client.computations.get_result_file_urls(computation_id, path).download_url
    with urllib.request.urlopen(url) as response:
        return response.read().decode()


def check_results(client: CodeOcean, computation_id: str, giant_sha: str) -> tuple[list[str], dict]:
    """Return failed-check messages and the parsed smoke-test summary."""
    failures = []
    sizes = {item.path: item.size for item in client.computations.list_computation_results(computation_id).items}

    for name in (COMMIT_FILE, TRIAL_TABLE_FILE, SUMMARY_FILE):
        if not sizes.get(name):
            failures.append(f"missing or empty result: {name}")

    if sizes.get(COMMIT_FILE):
        ran = read_result_text(client, computation_id, COMMIT_FILE).strip()
        if ran != giant_sha:
            failures.append(f"capsule ran GIAnT-MATLAB {ran}, expected {giant_sha}")

    summary = {}
    if sizes.get(SUMMARY_FILE):
        summary = json.loads(read_result_text(client, computation_id, SUMMARY_FILE))
    return failures, summary


def write_step_summary(lines: list[str]) -> None:
    """Append lines to the GitHub Actions job summary when running in Actions."""
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as f:
            f.write("\n".join(lines) + "\n")


def main() -> int:
    """Run the capsule, wait for it, and return a process exit code."""
    args = parse_args()
    client = CodeOcean(domain=CO_DOMAIN, token=os.environ["CODEOCEAN_TOKEN"])

    data_assets = None
    if args.session_asset_id:
        mount = client.data_assets.get_data_asset(args.session_asset_id).mount
        data_assets = [DataAssetsRunParam(id=args.session_asset_id, mount=mount)]
    computation = client.computations.run_capsule(
        RunParams(
            capsule_id=args.capsule_id,
            data_assets=data_assets,
            parameters=["--giant_sha", args.giant_sha],
        )
    )
    run_url = f"{CO_DOMAIN}/capsule/{client.capsules.get_capsule(args.capsule_id).slug}/tree"
    print(f"Started computation {computation.id} for GIAnT-MATLAB {args.giant_sha}", flush=True)

    try:
        computation = client.computations.wait_until_completed(
            computation, polling_interval=POLL_SECONDS, timeout=args.timeout_minutes * 60
        )
    except TimeoutError:
        failures = [f"computation did not finish within {args.timeout_minutes} minutes"]
        summary = {}
    else:
        failures = []
        if computation.state != ComputationState.Completed:
            failures.append(f"computation state is {computation.state.value}")
        if computation.end_status != ComputationEndStatus.Succeeded:
            status = computation.end_status.value if computation.end_status else None
            failures.append(f"computation end status is {status}")
        if computation.has_results:
            result_failures, summary = check_results(client, computation.id, args.giant_sha)
            failures.extend(result_failures)
        else:
            failures.append("computation produced no results")
            summary = {}

    verdict = "failed" if failures else "passed"
    lines = [
        f"## Code Ocean integration check {verdict}",
        f"- GIAnT-MATLAB commit: `{args.giant_sha}`",
        f"- Computation: `{computation.id}` ([capsule]({run_url}))",
    ]
    lines += [f"- Failure: {message}" for message in failures]
    if summary:
        lines += ["", "```json", json.dumps(summary, indent=2), "```"]
    write_step_summary(lines)
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

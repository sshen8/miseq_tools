import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

import requests

from .samplesheet import format_samplesheet
from .demux_stats import demux

DEFAULT_PARENT = "/n/liau_lab/NextGenRunData"
SBATCH_ARGS = ["-p", "test", "-c", "1", "-t", "12:00:00", "--mem", "32G", "-J", "miseq_download"]
RUN_FAILED = {"Failed", "Stopped", "Timed Out", "Needs Attention"}
APPSESSION_FAILED = {"Failed", "Aborted"}
MAX_BS_ERRORS = 3


# ---------------------------------------------------------------- helpers

def _bs(*args):
    """Run a `bs` command and return its parsed JSON output."""
    res = subprocess.run(["bs", *args, "-f", "json"], check=True, capture_output=True, text=True)
    return json.loads(res.stdout) if res.stdout.strip() else None


def _submit(state_path, delay=True):
    """Submit a job that runs one `download-step` on the state file, 5 minutes from now if `delay`."""
    state_path = Path(state_path).resolve()
    state = json.loads(state_path.read_text())
    wrap = f"{sys.executable} -m miseq_tools download-step {state_path}"
    begin = ["--begin=now+5minutes"] if delay else []
    res = subprocess.run(["sbatch", "--parsable", *SBATCH_ARGS, *begin, "-o", str(state_path.parent / "slurm-%j.out"), f"--wrap={wrap}"],
                         check=True, capture_output=True, text=True)
    job_id = res.stdout.strip().split(";")[0]
    state["jobs"].append(job_id)
    state_path.write_text(json.dumps(state, indent=2))
    logging.info(f"Submitted job {job_id}")
    return job_id

def _slack(webhook, text, blocks=None):
    """Post to Slack; `text` is the notification fallback when `blocks` are given."""
    payload = {"text": text}
    if blocks:
        payload["blocks"] = blocks
    if not webhook:
        logging.info(f"No Slack webhook, not sending: {json.dumps(payload)}")
        return
    try:
        requests.post(webhook, json=payload, timeout=30).raise_for_status()
    except requests.RequestException:
        logging.exception("Failed to send Slack message")


def _is_bclconvert(appsession):
    return (appsession.get("Application", {}).get("AppFamilySlug") == "illumina-inc.bcl-convert"
            or re.match(r"^BCL ?Convert", appsession.get("Name", ""), flags=re.IGNORECASE) is not None)


def _find_bclconvert(run_id):
    """Newest BCL Convert appsession for the run, or None."""
    appsessions = [a for a in (_bs("list", "appsessions", f"--input-run={run_id}") or []) if _is_bclconvert(a)]
    return max(appsessions, key=lambda a: a["DateCreated"], default=None)


def _find_stats(outdir):
    for fname in ("Demultiplex_Stats.csv", "Stats.json"):
        if found := sorted(Path(outdir).rglob(fname)):
            return found[0]
    return None


def _table_block(df_pool):
    """Slack Block Kit table block (max 100 rows, including the header)."""
    def bold(text):
        return {"type": "rich_text", "elements": [{"type": "rich_text_section", "elements": [
            {"type": "text", "text": text, "style": {"bold": True}}]}]}

    rows = [[bold("Pool label"), bold("Requested reads"), bold("Actual reads")]]
    for label, row in df_pool.iterrows():
        rows.append([{"type": "raw_text", "text": str(label)},
                     {"type": "raw_text", "text": f"{row['intended']:,.0f}"},
                     {"type": "raw_text", "text": f"{row['actual']:,.0f}"}])
    if len(rows) > 100:
        logging.warning(f"Slack table truncated to 99 of {len(rows) - 1} pools")
        rows = rows[:100]
    return {"type": "table", "column_settings": [{"is_wrapped": True}, {"align": "right"}, {"align": "right"}], "rows": rows}


# ---------------------------------------------------------------- slurm step

def step(state_path):
    """Advance the pipeline as far as possible; resubmit itself if BaseSpace isn't ready yet."""
    state_path = Path(state_path).resolve()
    state = json.loads(state_path.read_text())

    def save():
        state_path.write_text(json.dumps(state, indent=2))

    def fail(msg):
        logging.error(msg)
        state["phase"] = "failed"
        save()
        _slack(state.get("slack_webhook"), f"Run {state['run_name']} failed: {msg}\nLog: {state['outdir']}")

    def resubmit():
        save()
        _submit(state_path)

    if state["phase"] in ("done", "failed"):
        logging.info(f"Nothing to do, phase is {state['phase']}")
        return

    try:
        if state["phase"] == "wait_run":
            try:
                status = _bs("get", "run", "-i", state["run_id"])["Status"]
            except subprocess.CalledProcessError as e:
                return _bs_error(state, e, resubmit, fail)
            state["bs_errors"] = 0
            logging.info(f"Run {state['run_name']} status: {status}")
            if status in RUN_FAILED:
                return fail(f"run status is {status}")
            if status != "Complete":
                return resubmit()
            state["phase"] = "wait_appsession"
            save()

        if state["phase"] == "wait_appsession":
            try:
                appsession = _find_bclconvert(state["run_id"])
            except subprocess.CalledProcessError as e:
                return _bs_error(state, e, resubmit, fail)
            state["bs_errors"] = 0
            if appsession is None:
                logging.info("No BCL Convert appsession yet")
                return resubmit()
            status = appsession["ExecutionStatus"]
            logging.info(f"Appsession {appsession['Name']} ({appsession['Id']}) status: {status}")
            if status in APPSESSION_FAILED:
                return fail(f"BCL Convert appsession {appsession['Id']} status is {status}")
            if status != "Complete":
                return resubmit()
            state["appsession_id"] = appsession["Id"]
            state["phase"] = "download"
            save()

        if state["phase"] == "download":
            logging.info(f"Downloading appsession {state['appsession_id']} to {state['outdir']}")
            subprocess.run(["bs", "download", "appsession", "--id", state["appsession_id"], "-o", state["outdir"], "--no-progress-bars"], check=True)
            state["phase"] = "report"
            save()

        if state["phase"] == "report":
            text = f"Sequencing run is complete! Data is available at `{state['outdir']}`"
            blocks = None
            if state["samplesheet"]:
                stats = _find_stats(state["outdir"])
                if stats is None:
                    raise FileNotFoundError(f"No Demultiplex_Stats.csv or Stats.json in {state['outdir']}")
                df_pool = demux(state["samplesheet"], str(stats), out=str(Path(state["outdir"]) / "demux_stats.pdf"))
                blocks = [{"type": "section", "text": {"type": "mrkdwn", "text": text}}, _table_block(df_pool)]
            _slack(state.get("slack_webhook"), text, blocks)
            state["phase"] = "done"
            save()
    except Exception:
        traceback.print_exc()
        fail(f"unexpected error in phase {state['phase']}, see logs/slurm-*.out")


def _bs_error(state, e, resubmit, fail):
    state["bs_errors"] = state.get("bs_errors", 0) + 1
    logging.warning(f"bs failed ({state['bs_errors']}/{MAX_BS_ERRORS}): {e.stderr}")
    if state["bs_errors"] >= MAX_BS_ERRORS:
        return fail(f"bs failed {MAX_BS_ERRORS} times in a row: {e.stderr}")
    return resubmit()


# ---------------------------------------------------------------- TUI

def _build_app():
    from textual import work
    from textual.app import App, ComposeResult
    from textual.containers import Horizontal, VerticalScroll
    from textual.widgets import Button, Footer, Header, Input, Label, OptionList, RichLog
    from textual.widgets.option_list import Option

    class LogToWidget(logging.Handler):
        def __init__(self, widget):
            super().__init__()
            self.widget = widget

        def emit(self, record):
            color = {"WARNING": "yellow", "ERROR": "red", "CRITICAL": "red"}.get(record.levelname, "white")
            self.widget.write(f"[{color}]{record.levelname}: {record.getMessage()}[/]")

    class DownloadApp(App):
        TITLE = "MiSeq download"
        CSS = """
        .step { margin: 1 0 0 0; text-style: bold; }
        Horizontal { height: auto; }
        Input { width: 1fr; }
        OptionList { height: 12; }
        RichLog { height: 12; border: round $primary; }
        """
        BINDINGS = [("q", "quit", "Quit")]

        def __init__(self):
            super().__init__()
            self.tmpdir = tempfile.TemporaryDirectory()
            self.samplesheet = None  # validated input path, or None if skipped
            self.sheet_done = False
            self.outdir = None

        def compose(self) -> ComposeResult:
            yield Header()
            with VerticalScroll():
                yield Label("1. Sample sheet (optional)", classes="step")
                with Horizontal():
                    yield Input(placeholder="path/to/samplesheet.csv", id="sheet")
                    yield Button("Check", id="check")
                    yield Button("Skip", id="skip")
                yield Label("2. Output directory", classes="step")
                with Horizontal():
                    yield Input(f"{DEFAULT_PARENT}/{datetime.date.today():%y%m%d}_miseq", id="outdir", disabled=True)
                    yield Button("Create", id="create", disabled=True)
                yield Label("3. BaseSpace run", classes="step")
                yield OptionList(id="runs", disabled=True)
                with Horizontal():
                    yield Button("Refresh runs", id="refresh", disabled=True)
                yield Label("4. Slack webhook (optional; leave blank for no message)", classes="step")
                with Horizontal():
                    yield Input(os.environ.get("SLACK_WEBHOOK_URL"), placeholder="https://hooks.slack.com/services/...", id="webhook")
                    yield Button("Submit", id="submit", variant="success", disabled=True)
                yield RichLog(markup=True, wrap=True, id="log")
            yield Footer()

        def on_mount(self):
            self.log_widget = self.query_one("#log", RichLog)
            logging.getLogger().handlers.clear()  # stderr output would corrupt the TUI
            self.log_handler = LogToWidget(self.log_widget)
            logging.getLogger().addHandler(self.log_handler)
            logging.getLogger().setLevel(logging.INFO)
            if shutil.which("bs") is None:
                logging.error("BaseSpace CLI `bs` not found on PATH")
            if shutil.which("sbatch") is None:
                logging.error("`sbatch` not found on PATH")
            if not self.query_one("#webhook", Input).value:
                logging.warning("No SLACK_WEBHOOK_URL found in environment or .env; enter one in step 4 to get a Slack message")

        def on_unmount(self):
            logging.getLogger().removeHandler(self.log_handler)

        def _enable(self, *ids):
            for i in ids:
                self.query_one(f"#{i}").disabled = False

        def _sheet_done(self):
            self.sheet_done = True
            self._enable("outdir", "create")
            self.query_one("#outdir", Input).focus()

        def on_button_pressed(self, event: Button.Pressed):
            getattr(self, f"_on_{event.button.id}")()

        def _on_check(self):
            path = self.query_one("#sheet", Input).value.strip()
            if not path:
                return self._on_skip()
            path = Path(path).expanduser().resolve()
            logging.info(f"Checking {path}")
            try:
                format_samplesheet(str(path), str(Path(self.tmpdir.name) / "samplesheet.csv"))
            except Exception as e:
                logging.error(f"Could not read sample sheet: {e!r}")
                return
            logging.info("Sample sheet checked. Fix and re-check, or continue to step 2.")
            self.samplesheet = path
            self._sheet_done()

        def _on_skip(self):
            self.samplesheet = None
            logging.info("No sample sheet; demux stats will be skipped")
            self._sheet_done()

        def _on_create(self):
            outdir = Path(self.query_one("#outdir", Input).value.strip()).expanduser().resolve()
            try:
                outdir.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                logging.error(f"Could not create {outdir}: {e}")
                return
            if any(p for p in outdir.iterdir()):
                logging.warning(f"{outdir} is not empty")
            if self.samplesheet is not None:
                shutil.copy(self.samplesheet, outdir / "samplesheet_original.csv")
                shutil.copy(Path(self.tmpdir.name) / "samplesheet.csv", outdir / "samplesheet.csv")
            self.outdir = outdir
            logging.info(f"Using output directory {outdir}")
            self._enable("runs", "refresh")
            self._on_refresh()

        def _on_refresh(self):
            self.query_one("#runs", OptionList).clear_options()
            self.load_runs()

        @work(thread=True, exclusive=True)
        def load_runs(self):
            try:
                runs = sorted(_bs("list", "runs") or [], key=lambda r: r["DateCreated"], reverse=True)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                self.call_from_thread(logging.error, f"bs list runs failed: {getattr(e, 'stderr', e)}")
                return
            options = [Option(f"{r['Name']}  —  {r['Status']}  —  {r['DateCreated']}", id=r["Id"]) for r in runs]
            self.call_from_thread(self._show_runs, options)

        def _show_runs(self, options):
            runs = self.query_one("#runs", OptionList)
            runs.add_options(options)
            runs.focus()

        def on_option_list_option_highlighted(self, event: OptionList.OptionHighlighted):
            self._enable("submit")

        def _on_submit(self):
            runs = self.query_one("#runs", OptionList)
            if runs.highlighted is None:
                return
            option = runs.get_option_at_index(runs.highlighted)
            state = {
                "run_id": option.id,
                "run_name": str(option.prompt).split("  —  ")[0],
                "outdir": str(self.outdir),
                "samplesheet": str(self.outdir / "samplesheet_original.csv") if self.samplesheet else None,
                "slack_webhook": self.query_one("#webhook", Input).value.strip(),
                "phase": "wait_run",
                "jobs": [],
            }
            log_dir = self.outdir / "logs"
            log_dir.mkdir(exist_ok=True)  # slurm won't create it
            state_path = log_dir / "download_state.json"
            state_path.write_text(json.dumps(state, indent=2))
            try:
                job_id = _submit(state_path, delay=False)
            except subprocess.CalledProcessError as e:
                logging.error(f"sbatch failed: {e.stderr}")
                return
            self.exit(f"Submitted job {job_id} for run {state['run_name']}.\n"
                      f"It will resubmit itself every 5 minutes until the data is downloaded.\n"
                      f"State: {state_path}\nLogs: {log_dir}/slurm-*.out\n"
                      f"To stop: scancel the latest job in {state_path} (or `scancel -n miseq_download`).")

    return DownloadApp()


def run_tui():
    message = _build_app().run()
    if message:
        print(message)

import csv
import datetime
import json
import os
import threading
import time
from typing import Any, Dict, List, Optional


_try_import_psutil = True
try:
    import psutil
except ImportError:  # pragma: no cover
    psutil = None
    _try_import_psutil = False


class ResourceMonitor:
    """
    Track memory and CPU consumption for a block of pipeline execution.

    Usage:
        monitor = ResourceMonitor(interval=1.0, report_path="monitor_report.txt")
        monitor.start()
        run_pipeline()
        monitor.stop()
        monitor.write_report()

    Or as a context manager:
        with ResourceMonitor(interval=1.0, report_path="monitor_report.txt") as monitor:
            run_pipeline()
            # report is written automatically on exit if report_path is provided
    """

    def __init__(self, interval: float = 1.0, report_path: Optional[str] = None, label: Optional[str] = None):
        self.interval = float(interval)
        self.report_path = report_path
        self.label = label or "pipeline"
        self.samples: List[Dict[str, Any]] = []
        self._thread: Optional[threading.Thread] = None
        self._stop_event = threading.Event()
        self._started = False
        self._start_time: Optional[float] = None
        self._end_time: Optional[float] = None

    def _sample(self) -> None:
        timestamp = datetime.datetime.now().isoformat()
        cpu_percent = None
        memory_rss = None
        memory_vms = None
        memory_percent = None

        if psutil is not None:
            proc = psutil.Process(os.getpid())
            cpu_percent = proc.cpu_percent(interval=None)
            mem_info = proc.memory_info()
            memory_rss = mem_info.rss
            memory_vms = mem_info.vms
            memory_percent = proc.memory_percent()
        else:
            cpu_percent = self._sample_cpu_fallback()
            memory_rss, memory_vms, memory_percent = self._sample_memory_fallback()

        self.samples.append(
            {
                "timestamp": timestamp,
                "cpu_percent": cpu_percent,
                "memory_rss": memory_rss,
                "memory_vms": memory_vms,
                "memory_percent": memory_percent,
            }
        )

    def _sample_cpu_fallback(self) -> float:
        # Fallback: return 0.0 if we cannot sample reliably
        return 0.0

    def _sample_memory_fallback(self) -> tuple[Optional[float], Optional[float], Optional[float]]:
        rss = None
        vms = None
        mem_percent = None
        if os.path.exists("/proc/self/statm"):
            try:
                with open("/proc/self/statm", "r", encoding="utf-8") as handle:
                    parts = handle.read().split()
                if len(parts) >= 2:
                    page_size = os.sysconf("SC_PAGE_SIZE")
                    rss = int(parts[1]) * page_size
                    vms = int(parts[0]) * page_size
            except Exception:
                rss = None
                vms = None
        return rss, vms, mem_percent

    def _run(self) -> None:
        self._started = True
        self._start_time = time.time()
        if psutil is not None:
            psutil.Process(os.getpid()).cpu_percent(interval=None)
        while not self._stop_event.wait(self.interval):
            self._sample()
        self._sample()
        self._end_time = time.time()

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop_event.clear()
        self.samples = []
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop_event.set()
        self._thread.join(timeout=self.interval * 2 + 1)
        if self._thread.is_alive():
            raise RuntimeError("ResourceMonitor thread did not stop")
        if self.report_path:
            self.write_report(self.report_path)

    def __enter__(self) -> "ResourceMonitor":
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.stop()

    def _convert_bytes(self, value: Optional[float]) -> Optional[str]:
        if value is None:
            return None
        units = ["B", "KiB", "MiB", "GiB", "TiB"]
        amount = float(value)
        for unit in units:
            if amount < 1024.0 or unit == units[-1]:
                return f"{amount:.2f} {unit}"
            amount /= 1024.0
        return f"{amount:.2f} TiB"

    def get_summary(self) -> Dict[str, Any]:
        if not self.samples:
            return {
                "label": self.label,
                "duration_seconds": 0.0,
                "sample_count": 0,
                "cpu_avg_percent": 0.0,
                "cpu_peak_percent": 0.0,
                "memory_rss_avg_bytes": 0,
                "memory_rss_peak_bytes": 0,
                "memory_vms_avg_bytes": 0,
                "memory_vms_peak_bytes": 0,
                "memory_percent_avg": 0.0,
                "memory_percent_peak": 0.0,
            }

        cpu_values = [sample["cpu_percent"] for sample in self.samples if sample["cpu_percent"] is not None]
        mem_rss_values = [sample["memory_rss"] for sample in self.samples if sample["memory_rss"] is not None]
        mem_vms_values = [sample["memory_vms"] for sample in self.samples if sample["memory_vms"] is not None]
        mem_pct_values = [sample["memory_percent"] for sample in self.samples if sample["memory_percent"] is not None]

        duration = (self._end_time or time.time()) - (self._start_time or time.time())
        summary = {
            "label": self.label,
            "duration_seconds": round(duration, 2),
            "sample_count": len(self.samples),
            "cpu_avg_percent": round(sum(cpu_values) / len(cpu_values), 2) if cpu_values else 0.0,
            "cpu_peak_percent": round(max(cpu_values), 2) if cpu_values else 0.0,
            "memory_rss_avg_bytes": int(sum(mem_rss_values) / len(mem_rss_values)) if mem_rss_values else 0,
            "memory_rss_peak_bytes": int(max(mem_rss_values)) if mem_rss_values else 0,
            "memory_vms_avg_bytes": int(sum(mem_vms_values) / len(mem_vms_values)) if mem_vms_values else 0,
            "memory_vms_peak_bytes": int(max(mem_vms_values)) if mem_vms_values else 0,
            "memory_percent_avg": round(sum(mem_pct_values) / len(mem_pct_values), 2) if mem_pct_values else 0.0,
            "memory_percent_peak": round(max(mem_pct_values), 2) if mem_pct_values else 0.0,
        }
        summary.update(
            {
                "memory_rss_avg": self._convert_bytes(summary["memory_rss_avg_bytes"]),
                "memory_rss_peak": self._convert_bytes(summary["memory_rss_peak_bytes"]),
                "memory_vms_avg": self._convert_bytes(summary["memory_vms_avg_bytes"]),
                "memory_vms_peak": self._convert_bytes(summary["memory_vms_peak_bytes"]),
            }
        )
        return summary

    def get_report_dict(self) -> Dict[str, Any]:
        report = {
            "label": self.label,
            "started_at": datetime.datetime.now().isoformat(),
            "summary": self.get_summary(),
            "samples": self.samples,
        }
        return report

    def write_report(self, path: Optional[str] = None, fmt: str = "txt") -> None:
        if path is None:
            path = self.report_path
        if path is None:
            raise ValueError("A report path is required to write the report.")

        fmt = fmt.lower()
        if fmt not in {"txt", "csv", "json"}:
            raise ValueError("Report format must be one of: txt, csv, json")

        report = self.get_report_dict()
        directory = os.path.dirname(path)
        if directory and not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        if fmt == "json" or path.endswith(".json"):
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(report, handle, indent=2)
        elif fmt == "csv" or path.endswith(".csv"):
            self._write_csv(path)
        else:
            self._write_txt(path)

    def _write_txt(self, path: str) -> None:
        summary = self.get_summary()
        lines = [
            f"Resource report for: {summary['label']}",
            f"Duration (s): {summary['duration_seconds']}",
            f"Samples: {summary['sample_count']}",
            f"CPU average (%): {summary['cpu_avg_percent']}",
            f"CPU peak (%): {summary['cpu_peak_percent']}",
            f"Memory RSS average: {summary['memory_rss_avg']}",
            f"Memory RSS peak: {summary['memory_rss_peak']}",
            f"Memory VMS average: {summary['memory_vms_avg']}",
            f"Memory VMS peak: {summary['memory_vms_peak']}",
            f"Memory percent average: {summary['memory_percent_avg']}",
            f"Memory percent peak: {summary['memory_percent_peak']}",
            "",
            "Samples:",
        ]
        for sample in self.samples:
            lines.append(
                f"{sample['timestamp']}, cpu={sample['cpu_percent']}%, rss={self._convert_bytes(sample['memory_rss'])}, "
                f"vms={self._convert_bytes(sample['memory_vms'])}, percent={sample['memory_percent']}"
            )
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines))

    def _write_csv(self, path: str) -> None:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["timestamp", "cpu_percent", "memory_rss", "memory_vms", "memory_percent"],
            )
            writer.writeheader()
            for sample in self.samples:
                writer.writerow(sample)

    def print_summary(self) -> None:
        summary = self.get_summary()
        print("Resource consumption report:")
        print(f"  Label: {summary['label']}")
        print(f"  Duration: {summary['duration_seconds']} s")
        print(f"  Samples: {summary['sample_count']}")
        print(f"  CPU avg: {summary['cpu_avg_percent']}%")
        print(f"  CPU peak: {summary['cpu_peak_percent']}%")
        print(f"  Memory RSS avg: {summary['memory_rss_avg']}")
        print(f"  Memory RSS peak: {summary['memory_rss_peak']}")
        print(f"  Memory percent avg: {summary['memory_percent_avg']}%")
        print(f"  Memory percent peak: {summary['memory_percent_peak']}%")


def monitor_function(
    func,
    *args,
    report_path: Optional[str] = None,
    interval: float = 1.0,
    label: Optional[str] = None,
    **kwargs,
) -> tuple[Any, ResourceMonitor]:
    monitor = ResourceMonitor(interval=interval, report_path=report_path, label=label)
    monitor.start()
    try:
        result = func(*args, **kwargs)
    finally:
        monitor.stop()
    return result, monitor

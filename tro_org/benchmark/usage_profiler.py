import psutil
import os
import time
import threading
from datetime import datetime
import functools

USE_NVML = False
try:
    import pynvml
    pynvml.nvmlInit()
    USE_NVML = True
except Exception:
    USE_NVML = False

LOG_FILE = "resource_usage.tsv"

def get_vram_used_mb():
    """Return total VRAM used (sum over all GPUs) in MB. 0.0 if NVML not available."""
    if not USE_NVML:
        return 0.0

    try:
        device_count = pynvml.nvmlDeviceGetCount()
        total_used = 0
        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
            total_used += mem.used  # bytes
        return total_used / (1024 ** 2)
    except Exception:
        return 0.0


def get_gpu_util_percent():
    """
    Return total GPU utilization (%) summed over all GPUs.
    On systems without NVIDIA+NVML, returns 0.0.
    """
    if not USE_NVML:
        return 0.0

    try:
        device_count = pynvml.nvmlDeviceGetCount()
        total_util = 0
        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            util = pynvml.nvmlDeviceGetUtilizationRates(handle)
            total_util += util.gpu  # gpu utilization in %
        return float(total_util)
    except Exception:
        return 0.0

def get_gpu_memory_mb():
    """
    Return total GPU memory used (sum over all GPUs) in MB.
    Returns 0.0 if NVML is unavailable (e.g., on M1 Mac).
    """
    if not USE_NVML:
        return 0.0

    try:
        total_used = 0
        device_count = pynvml.nvmlDeviceGetCount()
        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            mem = pynvml.nvmlDeviceGetMemoryInfo(handle)
            total_used += mem.used  # bytes
        return total_used / (1024 ** 2)
    except:
        return 0.0



def profile_resources(method_name: str, log_file: str = None, interval: float = 1.0):
    if log_file is None:
        log_file = LOG_FILE

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            process = psutil.Process(os.getpid())
            stop_event = threading.Event()

            os.makedirs(os.path.dirname(log_file), exist_ok=True)

            if not os.path.exists(LOG_FILE):
                with open(log_file, "w") as f:
                    f.write(
                        "method_name\t"
                        "timestamp\t"
                        "seconds_since_start\t"
                        "ram_mb\t"
                        "virt_mb\t"
                        "cpu_percent\t"
                        "gpu_util_percent\t"
                        "gpu_memory_mb\n"
                    )
            else:
                with open(log_file, "w") as f:
                    f.write(
                        "method_name\t"
                        "timestamp\t"
                        "seconds_since_start\t"
                        "ram_mb\t"
                        "virt_mb\t"
                        "cpu_percent\t"
                        "gpu_util_percent\t"
                        "gpu_memory_mb\n"
                    )
                    f.close()

            def monitor():
                start = time.time()
                while not stop_event.is_set():
                    now = time.time()
                    elapsed = now - start

                    mem_info = process.memory_info()
                    ram_mb = mem_info.rss / (1024 ** 2)
                    virt_mb = mem_info.vms / (1024 ** 2)
                    cpu_percent = process.cpu_percent(interval=None)
                    gpu_util = get_gpu_util_percent()
                    gpu_mem_mb = get_gpu_memory_mb()

                    with open(log_file, "a") as f:
                        f.write(
                            f"{method_name}\t"
                            f"{datetime.now().isoformat()}\t"
                            f"{elapsed:.1f}\t"
                            f"{ram_mb:.2f}\t"
                            f"{virt_mb:.2f}\t"
                            f"{cpu_percent:.1f}\t"
                            f"{gpu_util:.1f}\t"
                            f"{gpu_mem_mb:.2f}\n"
                        )

                    stop_event.wait(interval)

            t = threading.Thread(target=monitor, daemon=True)
            t.start()
            try:
                return func(*args, **kwargs)
            finally:
                stop_event.set()
                t.join(timeout=2.0)

        return wrapper
    return decorator
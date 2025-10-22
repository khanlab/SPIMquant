import sys
import threading
import traceback

def dump_all_thread_stacks():
    for thread_id, frame in sys._current_frames().items():
        thread = next((t for t in threading.enumerate() if t.ident == thread_id), None)
        name = thread.name if thread else "unknown"
        print(f"\n--- Thread: {name} ({thread_id}) ---")
        traceback.print_stack(frame)


import signal

def handler(signum, frame):
    print("Received SIGUSR1, dumping thread stacks...")
    dump_all_thread_stacks()

signal.signal(signal.SIGUSR1, handler)


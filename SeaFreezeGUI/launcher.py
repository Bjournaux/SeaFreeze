"""Desktop launcher for SeaFreeze GUI.

Runs Streamlit in-process — no subprocess, no external Python needed.
The PyInstaller bundle includes the full Python interpreter + all deps.
"""

import os
import sys
import socket
import threading
import time
import webbrowser


def _resource_path(relative):
    """Resolve path whether running from source or PyInstaller bundle."""
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, relative)
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), relative)


def _find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _open_browser_when_ready(port, timeout=30):
    """Background thread: wait for the server, then open the browser once."""
    url = f"http://localhost:{port}"
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                webbrowser.open(url)
                print(f"SeaFreeze GUI running at {url}")
                return
        except OSError:
            time.sleep(0.5)
    print("Warning: timed out waiting for Streamlit to start.", file=sys.stderr)


def main():
    port = _find_free_port()
    app_path = _resource_path("app.py")

    # Point Streamlit at the bundled .streamlit config
    config_dir = _resource_path(".streamlit")
    if os.path.isdir(config_dir):
        os.environ["STREAMLIT_CONFIG_DIR"] = config_dir

    # Streamlit config via env vars
    os.environ["STREAMLIT_GLOBAL_DEVELOPMENT_MODE"] = "false"
    os.environ["STREAMLIT_SERVER_PORT"] = str(port)
    os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"
    os.environ["STREAMLIT_THEME_PRIMARY_COLOR"] = "#1f77b4"

    # Open browser in a background thread once the server is ready
    browser_thread = threading.Thread(
        target=_open_browser_when_ready,
        args=(port,),
        daemon=True,
    )
    browser_thread.start()

    # Run Streamlit in-process (no subprocess needed)
    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit", "run", app_path,
        "--server.port", str(port),
        "--server.headless", "true",
    ]

    print(f"Starting SeaFreeze GUI on port {port}...")
    print("Close this window or press Ctrl+C to stop.\n")

    stcli.main()


if __name__ == "__main__":
    main()

"""Entry point: `plviewer` launches the Streamlit app."""
import subprocess
import sys
from pathlib import Path

APP_PATH = Path(__file__).resolve().parent.parent / "app.py"


def main() -> None:
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(APP_PATH), *sys.argv[1:]], check=True)


if __name__ == "__main__":
    main()

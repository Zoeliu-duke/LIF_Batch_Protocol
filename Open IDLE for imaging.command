#!/bin/zsh
# Double-click in Finder: opens IDLE with the imaging Python (numpy, readlif, ... installed).
# Then File > Open... your process_lif.py and Run > Run Module (fn+F5).
if [ ! -x ~/.venvs/imaging/bin/python ]; then
  echo "The imaging Python (~/.venvs/imaging) is missing. See README.md, section 3 (One-time setup)."
  read -k1 "?Press any key to close."
  exit 1
fi
nohup ~/.venvs/imaging/bin/python -m idlelib >/dev/null 2>&1 &
sleep 1
osascript -e 'tell application "Terminal" to close (every window whose name contains "Open IDLE for imaging")' >/dev/null 2>&1 &
exit 0

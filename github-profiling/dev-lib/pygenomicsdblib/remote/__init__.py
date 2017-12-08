import sys
import os.path
import platform

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not root_path in sys.path:
    sys.path.append(root_path)

hostname = platform.node().split('.')[0]
__all__ = ["exec_test", "prerun_check", "stat_capturer"]

import sys
import os.path

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not root_path in sys.path:
    sys.path.append(root_path)
import shared.profiling_env as ppenv
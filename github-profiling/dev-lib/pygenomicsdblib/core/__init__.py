import sys
import os.path

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not root_path in sys.path:
    sys.path.append(root_path)

__all__ = ['loader_profiler']

from os.path import dirname, abspath, join, normpath
import sys

current_dir = dirname(abspath(__file__))
data_dir = join(current_dir,'../data')
sys.path.append(data_dir)

__all__ = ["config", "pathway", "reaction", "drawpathway"]
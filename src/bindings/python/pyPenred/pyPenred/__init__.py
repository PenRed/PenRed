# pyPenred/__init__.py
import importlib

simulator = importlib.import_module(".simulator", __name__)

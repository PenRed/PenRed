# Installing pyPenred

Certain features, such as DICOM loading or running simulations, require the `pyPenred` package inside Blender's Python environment. This section explains how to install and manage the package. Note that Blender uses its own embedded Python environment, which **is separate from your system's default Python installation**.

## Automatic Installation

If `pyPenred` is not installed, a warning banner will appear in any panel that depends on it:

<img src="../images/dependencyWarning.png" alt="pyPenred missing warning" width="500" style="display: block; margin: 0 auto"/>

To install the package automatically:

1. Click the **Open pyPenred Preferences** button in the warning banner to open the preferences window.
2. Click the **Install pyPenred** button:

<img src="../images/installPyPenredAddon.png" alt="pyPenred preferences panel" width="500" style="display: block; margin: 0 auto"/>

From this same window, you can also update or remove `pyPenred` at any time.

## Manual / Custom Installation

To install a specific or local version of `pyPenred`, you must use Blender's embedded Python executable.

1. Open Blender's **Scripting** workspace.
2. Run the following code in the Python Console:

```python
import sys
print(sys.executable)
```

3. Open your system terminal (Command Prompt/PowerShell on Windows, Terminal on Linux/macOS) and run pip using the exact path returned in step 2:

```bash
# Example (Linux):
/path/to/blender/python -m pip install pyPenred

# Example (Windows):
"C:\Path\To\Blender\python.exe" -m pip install pyPenred
```

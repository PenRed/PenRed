# Installing pyPenred

Certain features, such as DICOM loading or running simulations, require the `pyPenred` package inside Blender's Python environment. This section explains how to manage the package. Note that Blender uses its own embedded Python environment, which **is separate from your system's default Python installation**.

## Automatic Installation

The `pyPenred` module and all its dependencies are already bundled within the Blender plugin ZIP file, so no additional steps are required to enable it. 

Each plugin release bundles a specific `pyPenred` version. To use a different version or update it, check the **[PenRed releases page](https://github.com/PenRed/PenRed/releases)**.

If `pyPenred` is uninstalled or missing, a warning banner will appear in any panel that depends on it:

<img src="../images/dependencyWarning.png" alt="pyPenred missing warning" width="500" style="display: block; margin: 0 auto"/>

To resolve this issue, simply reinstall the Blender plugin.

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

**Note**: This method may fail due to a lack of write permissions, especially on macOS or system-wide Blender installations.

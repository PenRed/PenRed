import importlib
import importlib.util
import os
import shutil
import site
import subprocess
import sys

PACKAGE_NAME = "pyPenred"
MODULE_NAME = "pyPenred"  # Import name in Python


def get_user_site_packages():
    """Ensures user site-packages is in sys.path and returns the path."""
    user_site = site.getusersitepackages()
    if user_site not in sys.path:
        sys.path.append(user_site)
    return user_site


def is_installed():
    """Checks if the module can be found in sys.path without importing it."""
    get_user_site_packages()
    return importlib.util.find_spec(MODULE_NAME) is not None


def install_package():
    """Installs the package to the user site directory."""
    target_dir = get_user_site_packages()
    
    # pyYAML
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "pyyaml",
        "--target",
        target_dir,
        "--no-cache-dir"
    ]
    subprocess.check_call(cmd)

    # numpy
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "numpy",
        "--target",
        target_dir,
        "--no-cache-dir"
    ]
    subprocess.check_call(cmd)
    
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "-i",
        "https://test.pypi.org/simple/", # For test version
        "--extra-index-url",
        "https://pypi.org/simple",        
        PACKAGE_NAME,
        "--target",
        target_dir,
        "--no-cache-dir",
    ]
    subprocess.check_call(cmd)
    importlib.invalidate_caches()


def update_package():
    """Force-upgrades the package to the latest version."""
    target_dir = get_user_site_packages()
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--upgrade",
        PACKAGE_NAME,
        "--target",
        target_dir,
        "--no-cache-dir",
    ]
    subprocess.check_call(cmd)
    importlib.invalidate_caches()


def uninstall_package():
    """Uninstalls the package from the user environment."""
    target_dir = get_user_site_packages()

    # Call pip uninstall
    cmd = [sys.executable, "-m", "pip", "uninstall", "-y", PACKAGE_NAME]
    subprocess.call(cmd)

    # Fallback: remove files manually if installed inside target directory
    for item in os.listdir(target_dir):
        if item.lower().startswith(MODULE_NAME.lower()):
            item_path = os.path.join(target_dir, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
            elif os.path.isfile(item_path):
                os.remove(item_path)

    importlib.invalidate_caches()

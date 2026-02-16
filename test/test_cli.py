from pathlib import Path
import subprocess
import tomllib
import unittest

PROJECT_ROOT = Path(__file__).resolve().parent.parent 

def _get_binaries_from_makefile(makefile):
    for line in makefile.read_text().splitlines():
        if line.startswith("BINARIES"):
            return line.split("=")[1].split()
    raise ValueError("Failed to parse binaries from Makefile")

class TestCli(unittest.TestCase):
    """
    Minimal test case that each of the CLI entry points can be called from terminal
    """

    def test_cmd_gets_help(self):
        cfg = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
        for tool in cfg["project"]["scripts"]:
            result = subprocess.run(
                [tool, "-h"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertTrue(result.stdout.startswith(f"usage: {tool}"))

    def test_script_gets_help(self):
        scripts = [s.name for s in (PROJECT_ROOT / "scripts").iterdir()]
        for script in scripts:
            result = subprocess.run(
                [script, "-h"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertTrue(result.stderr.startswith(script))

    def test_binaries(self):
        binaries = _get_binaries_from_makefile(PROJECT_ROOT / "Makefile")
        for bn in binaries:
            result = subprocess.run(
                [bn, "--help"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if bn == "k8":
                # k8 writes it's help message to stderr not stdout
                # stdout has loads of javasript options
                self.assertTrue(bn in result.stderr)    
            else: 
                self.assertTrue(bn in result.stdout)
                self.assertFalse(result.stderr)


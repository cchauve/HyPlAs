import subprocess
import sys
import os

from pathlib import Path
from glob import glob
from typing import Optional, List, Dict

USE_VENV=True
DRY_RUN=False

class tool:
    name: str
    install: tuple[str, str, Optional[str]]
    depends: List["tool"]
    def __init__ (self, name, install, depends):
        self.name = name
        self.install= install
        self.depends = depends

    
class toolset:
    tools : Dict[str,"tool"]

    def __init__(self, tools):
        self.tools = {t[0] : tool(t[0], t[1], t[2]) for t in tools}

    def install(self):
        for name, t in self.tools.items():
            print(f"Installing {name}", file=sys.stderr)
            if t.install[0]:
                cmd = CMD(t.install[0], "install", t.install[1] if t.install[1] else name)
                ret = subprocess.run(cmd)
            else:
                t.install[2]()

def CMD(*args):
    cmd = []
    if DRY_RUN:
        cmd.append("echo")
    for a in args:
        cmd.append(a)
    return cmd
def install_prodigal_function():
    cmd = CMD(
        "git",
        "clone",
        "https://github.com/hyattpd/Prodigal.git"
        )
    subprocess.run(cmd)
    
    cwd = os.getcwd()
    os.chdir('Prodigal')

    cmd = CMD("git",
        "checkout",
        "GoogleImport"
        )
    
    subprocess.run(cmd)

    cmd = CMD("make")
    subprocess.run(cmd)
    
    cmd = CMD("make",
        "install",
        f"INSTALLDIR={ENV_DIR}/bin"
    )
    subprocess.run(cmd)
    
    os.chdir(cwd)

def install_spades_function():
    cmd = CMD("wget", "https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz")
    subprocess.run(cmd)
    cmd = CMD("tar", "-xzf", "SPAdes-4.2.0-Linux.tar.gz")
    subprocess.run(cmd)
    cmd = CMD("cp", "-rf")
    for file in glob("SPAdes-4.2.0-Linux/*"):
        cmd.append(file)
    cmd.append(f"{ENV_DIR}")
    subprocess.run(cmd)

def install_hyplass_function():
    cmd = CMD("make")
    subprocess.run(cmd)
    
    cmd = CMD("cp",
       "src/innotin",
       "src/select_missing_reads",
       "src/split_plasmid_reads",
       f"{ENV_DIR}/bin"
        )
    subprocess.run(cmd)

def install_minigraph_function():
    cmd = CMD("git", "clone", "https://github.com/lh3/minigraph.git")
    subprocess.run(cmd)

    cwd = os.getcwd()
    os.chdir('minigraph')

    cmd = CMD("git", "checkout", "v0.21")
    subprocess.run(cmd)

    cmd = CMD("make")
    subprocess.run(cmd)

    cmd = CMD("cp", "minigraph", f"{ENV_DIR}/bin")
    subprocess.run(cmd)

    os.chdir(cwd)

def install_racon_function():
    cmd = CMD("git", "clone", "--recursive", "https://github.com/lbcb-sci/racon.git", "racon")
    subprocess.run(cmd)

    cmd = CMD("mkdir", "-p", "racon/build")
    subprocess.run(cmd)

    cwd = os.getcwd()
    os.chdir("racon/build")
    
    cmd = CMD("cmake", "-DCMAKE_BUILD_TYPE=Release", "..")
    subprocess.run(cmd)

    cmd = CMD("make")
    subprocess.run(cmd)

    cmd = CMD("cp", "bin/racon", f"{ENV_DIR}/bin")
    subprocess.run(cmd)

if __name__ == "__main__":
    global ENV_DIR
    ENV_DIR = str(Path(sys.argv[1]).resolve())
    t = toolset([
        ["pandas", ("pip",None, None), []],
        ["numpy", ("pip",None,None), []],
        ["packaging", ("pip", None, None), []],
        ["platon", ("pip","cb-platon",None), ["prodigal"]],
        ["prodigal", (None, None, install_prodigal_function), []],
        ["spades", (None, None, install_spades_function), []], 
        ["racon", (None, None, install_racon_function), []],
        ["unicycler", ("pip", "git+https://github.com/rrwick/Unicycler.git", None), ["spades", "racon", "blast+"]],
        ["hyplass", (None, None, install_hyplass_function), ["unicycler", "platon", "packaging", "numpy", "pandas"]],
        ["minigraph", (None, None, install_minigraph_function), []]
    ])
    print("Nice")
    t.install()

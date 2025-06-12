import subprocess
import sys
import os
import argparse

from pathlib import Path
from glob import glob
from typing import Optional, List, Dict, Union
from types import FunctionType



def pip_install (*args):
    cmd = CMD("pip", "install")
    for a in args:
        cmd.append(a)
    return subprocess.run(cmd)

installers = {
    "pip": pip_install
}



class tool:
    name: str
    installation: Union[str, list, FunctionType]
    depends: List["tool"]
    done: bool
    def __init__ (self, name, installation, depends):
        self.name = name
        self.installation= installation
        self.depends = depends
        self.done = False
    def install(self, tools_dict):

        if self.done:
            return
        for tool in self.depends:
            tools_dict[tool].install(tools_dict)
        ins = self.installation

        print(f"\n#Installing {self.name}", file=sys.stderr)
        if not isinstance(ins, list):
            ins = [ins, self.name]
        if isinstance(ins[0], str):
            installers[ins[0]](*ins[1:])
        else:
            ins[0](*ins[1:])

        self.done = True
    
class toolset:
    tools : Dict[str,"tool"]

    def __init__(self, *tools):
        self.tools = {t[0] : tool(t[0], t[1], t[2]) for t in tools}

    def install(self, tool_name):
        assert tool_name in self.tools
        self.tools[tool_name].install(self.tools)

def change_dir(target_dir):
    if DRY_RUN:
        print(f"cd {target_dir}", file=sys.stderr)
    else:
        return os.chdir(target_dir)
def CMD(*args):
    cmd = []
    if DRY_RUN:
        cmd.append("echo")
    for a in args:
        cmd.append(a)
    return cmd
def install_prodigal_function(*x):
    cmd = CMD(
        "git",
        "clone",
        "https://github.com/hyattpd/Prodigal.git"
        )
    subprocess.run(cmd)
    
    cwd = os.getcwd()
    change_dir('Prodigal')

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
    
    change_dir(cwd)

def install_spades_function(*x):
    cmd = CMD("wget", "https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz")
    subprocess.run(cmd)
    cmd = CMD("tar", "-xzf", "SPAdes-4.2.0-Linux.tar.gz")
    subprocess.run(cmd)
    cmd = CMD("cp", "-rf")
    for file in glob("SPAdes-4.2.0-Linux/*"):
        cmd.append(file)
    cmd.append(f"{ENV_DIR}")
    subprocess.run(cmd)

def install_hyplass_function(*x):
    

    current_dir = os.getcwd()
    change_dir(SCRIPT_DIR)
    cmd = CMD("make")
    subprocess.run(cmd)
    


    cmd = CMD("cp",
       "src/innotin",
       "src/select_missing_reads",
       "src/split_plasmid_reads",
       f"{ENV_DIR}/bin"
        )
    subprocess.run(cmd)
    
    change_dir(current_dir)

def install_minigraph_function(*x):
    cmd = CMD("git", "clone", "https://github.com/lh3/minigraph.git")
    subprocess.run(cmd)

    cwd = os.getcwd()
    change_dir('minigraph')

    cmd = CMD("git", "checkout", "v0.21")
    subprocess.run(cmd)

    cmd = CMD("make")
    subprocess.run(cmd)

    cmd = CMD("cp", "minigraph", f"{ENV_DIR}/bin")
    subprocess.run(cmd)

    change_dir(cwd)

def install_racon_function(*x):
    cmd = CMD("git", "clone", "--recursive", "https://github.com/lbcb-sci/racon.git", "racon")
    subprocess.run(cmd)

    cmd = CMD("mkdir", "-p", "racon/build")
    subprocess.run(cmd)

    cwd = os.getcwd()
    change_dir("racon/build")
    
    cmd = CMD("cmake", "-DCMAKE_BUILD_TYPE=Release", "..")
    subprocess.run(cmd)

    cmd = CMD("make")
    subprocess.run(cmd)

    cmd = CMD("cp", "bin/racon", f"{ENV_DIR}/bin")
    subprocess.run(cmd)
    
    change_dir(cwd)


def installation_not_implemented(*args):
    print(f"Installation not implemented: {args}", file=sys.stderr)
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Hyplass Build Script")

    parser.add_argument("env", help="Path to virtual environment")
    parser.add_argument("-n", "--dry-run", help="Dry run", action="store_true")
    parser.add_argument("--tool", help="Tool to install", default="hyplass")
    parser.add_argument("--working-directory", help="Working directory for building tools", default=".")
    parser.add_argument("--skip", help="Installations to be skipped (i.e. prefer using available installations in the system", 
                        nargs="*", default=["spades", "blast+", "diamond", "mummer", "hmmer", "infernal", "minimap2"])
    args = parser.parse_args()

    global ENV_DIR
    global DRY_RUN
    global SCRIPT_DIR
    global TO_BE_SKIPPED
    DRY_RUN = args.dry_run
    ENV_DIR = str(Path(args.env).resolve())
    SCRIPT_DIR = os.getcwd()

    if args.working_directory != ".":
        subprocess.run(CMD("mkdir", "-p", args.working_directory))
        change_dir(args.working_directory)


    t = toolset(
        ["hyplass", install_hyplass_function, ["unicycler", "platon", "packaging", "numpy", "pandas", "minigraph", "minimap2"]],
        ["pandas", "pip", ["numpy"]],
        ["numpy", "pip", []],
        ["packaging", "pip", []],
        ["platon", ["pip","cb-platon"], ["prodigal", "diamond", "blast+", "mummer", "hmmer", "infernal"]],
        ["prodigal", install_prodigal_function, []],
        ["spades", install_spades_function, []], 
        ["racon", install_racon_function, []],
        ["unicycler", ["pip", "git+https://github.com/f0t1h/Unicycler.git"], ["spades", "racon"]],
        ["minigraph", install_minigraph_function, []],
        ["minimap2",installation_not_implemented, []],
        ["blast+",installation_not_implemented, []],
        ["diamond",installation_not_implemented, []],
        ["mummer",installation_not_implemented, []],
        ["hmmer",installation_not_implemented, []],
        ["infernal",installation_not_implemented, []],
    )

    for tool in args.skip:
        t.tools[tool].done = True

    t.install(args.tool)

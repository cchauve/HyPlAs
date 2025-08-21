import subprocess
import sys
import os
import argparse

from pathlib import Path
from glob import glob
from typing import Optional, List, Dict, Union
from types import FunctionType

def CMD(*args):
    cmd = []
    if DRY_RUN:
        cmd.append("echo")
    for a in args:
        cmd.append(a)
    return cmd

def run_or_exit(cmd, **kwargs):
    print(f"Running {' '.join(cmd)}", file=sys.stderr)
    comp_proc = subprocess.run(cmd)
    comp_proc.check_returncode()
    return comp_proc.returncode

def pip_install (*args):
    
    cmd = CMD("pip", "install")
    for a in args:
        cmd.append(a)
    return run_or_exit(cmd)

    

installers = {
    "pip": pip_install
}

def wget_link(link, *args):
    tgt = f"{link[link.rfind('/')+1:]}"
    print(tgt)
    if not os.path.isfile(tgt):
        cmd = CMD(
            "wget",
            link,
            *args
        )
        run_or_exit(cmd)
    return 0
def git_clone(repo, *args):
    tgt = f"{repo[repo.rfind('/')+1:repo.rfind('.git')]}"
    print(tgt)
    if not os.path.isdir(tgt):
        cmd = CMD(
            "git",
            "clone",
            repo,
            *args
        )
        run_or_exit(cmd)
    return 0

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
        print(f"cd {target_dir}", file=sys.stdout)
    else:
        return os.chdir(target_dir)

def install_prodigal_procedure(*x):
    git_clone(
        "https://github.com/hyattpd/Prodigal.git"
    )

    cwd = os.getcwd()
    change_dir('Prodigal')

    cmd = CMD("git",
        "checkout",
        "GoogleImport"
        )
    
    run_or_exit(cmd)

    cmd = CMD("make")
    run_or_exit(cmd)
    
    cmd = CMD("make",
        "install",
        f"INSTALLDIR={ENV_DIR}/bin"
    )
    run_or_exit(cmd)
    
    change_dir(cwd)


def install_blastn_procedure(*x):
    import platform
    if platform.system() == "Linux":
        os = "linux"
    elif platform.system() == "Darwin":
        os = "macosx"
    if(platform.machine() == 'arm64'):
        arch = "aarch64"
    else:
        if(os == "linux"):
            arch = "universal"
        else:
            arch = "x64"
    archive_file = f"ncbi-blast-2.17.0+-{arch}-{os}.tar.gz"
    
    cmd = wget_link(f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/{archive_file}")
   
    cmd = CMD("tar", "-xzf", archive_file)
    run_or_exit(cmd)
    cmd = CMD("cp", "-rf")
    for file in glob(f"{archive_file[:archive_file.rfind('tar.gz')]}/*"):
        cmd.append(file)
    cmd.append(f"{ENV_DIR}")
    run_or_exit(cmd)


def install_spades_procedure(*x):
    wget_link("https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz")
    cmd = CMD("tar", "-xzf", "SPAdes-4.2.0-Linux.tar.gz")
    run_or_exit(cmd)
    cmd = CMD("cp", "-rf")
    for file in glob("SPAdes-4.2.0-Linux/*"):
        cmd.append(file)
    cmd.append(f"{ENV_DIR}")
    run_or_exit(cmd)

def install_HyPlAs_cpp_binaries_procedure(*x):
    current_dir = os.getcwd()
    change_dir(SCRIPT_DIR)
    cmd = CMD("make")
    run_or_exit(cmd)
    
    cmd = CMD("cp",
       "src/innotin",
       "src/select_missing_reads",
       "src/split_plasmid_reads",
       f"{ENV_DIR}/bin"
        )
    run_or_exit(cmd)
        
    change_dir(current_dir)
def install_HyPlAs_procedure(*x):
    
    current_dir = os.getcwd()
    change_dir(SCRIPT_DIR)
    
    cmd = CMD("pip", "install", ".", "vvv", "--no-deps",  "--no-cache-dir")
    run_or_exit(cmd)
    
    change_dir(current_dir)

def install_minimap2_procedure(*x):
    import platform
    git_clone("https://github.com/lh3/minimap2.git")
    # run_or_exit(cmd)

    cwd = os.getcwd()
    change_dir('minimap2')
    
    cmd = CMD("git", "checkout", "v2.30")
    run_or_exit(cmd)

    if platform.machine() == "arm64":
        cmd = CMD("make", "arm_neon=1", "aarch64=1")
    else:
        cmd = CMD("make")
    run_or_exit(cmd)

    cmd = CMD("cp", "minimap2", f"{ENV_DIR}/bin")
    run_or_exit(cmd)

    change_dir(cwd)


def install_minigraph_procedure(*x):
    cmd = git_clone( "https://github.com/lh3/minigraph.git")

    cwd = os.getcwd()
    change_dir('minigraph')

    cmd = CMD("git", "checkout", "v0.21")
    run_or_exit(cmd)

    cmd = CMD("make")
    run_or_exit(cmd)

    cmd = CMD("cp", "minigraph", f"{ENV_DIR}/bin")
    run_or_exit(cmd)

    change_dir(cwd)

def install_racon_procedure(*x):
    cmd = git_clone("https://github.com/lbcb-sci/racon.git", "--recursive")

    cwd = os.getcwd()
    change_dir("racon")
    
    cmd = CMD("git", "checkout", "1.5.0")
    run_or_exit(cmd)

    cmd = CMD("mkdir", "-p", "build")
    run_or_exit(cmd)

    change_dir("build")
    

    cmd = CMD("cmake", "-DCMAKE_BUILD_TYPE=Release", "..")
    run_or_exit(cmd)

    cmd = CMD("make")
    run_or_exit(cmd)

    cmd = CMD("cp", "bin/racon", f"{ENV_DIR}/bin")
    run_or_exit(cmd)
    
    change_dir(cwd)


def installation_not_implemented(*args):
    print(f"Installation not implemented: {args}", file=sys.stderr)
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="HyPlAs Build Script")

    parser.add_argument("env", help="Path to virtual environment")
    parser.add_argument("-n", "--dry-run", help="Dry run", action="store_true")
    parser.add_argument("--tool", help="Tool to install", default="HyPlAs")
    parser.add_argument("--working-directory", help="Working directory for building tools", default="tmp-build")
    parser.add_argument("--skip", help="Installations to be skipped (i.e. prefer using available installations in the system", 
                        nargs="*", default=[  "diamond", "mummer", "hmmer", "infernal"])
    global args  
    args = parser.parse_args()
    global ENV_DIR
    global DRY_RUN
    global SCRIPT_DIR
 
    DRY_RUN = args.dry_run
    ENV_DIR = str(Path(args.env).resolve())
    SCRIPT_DIR = os.getcwd()

    if args.working_directory != ".":
        run_or_exit(CMD("mkdir", "-p", args.working_directory))
        change_dir(args.working_directory)


    t = toolset(
        ["HyPlAs", install_HyPlAs_procedure, ["HyPlAs_bin", "unicycler", "platon", "packaging", "numpy", "pandas", "minigraph", "minimap2"]],
        ["HyPlAs_bin", install_HyPlAs_cpp_binaries_procedure, []],
        ["pandas", "pip", ["numpy"]],
        ["numpy", "pip", []],
        ["packaging", "pip", []],
        ["platon", ["pip","cb-platon"], ["prodigal", "diamond", "mummer", "hmmer", "infernal"]],
        ["prodigal", install_prodigal_procedure, []],
        ["spades", install_spades_procedure, []], 
        ["racon", install_racon_procedure, []],
        ["unicycler", ["pip", "git+https://github.com/f0t1h/Unicycler.git"], ["spades", "racon"]],
        ["minigraph", install_minigraph_procedure, []],
        ["minimap2",install_minimap2_procedure, []],
        ["blast+",install_blastn_procedure, []],
        ["diamond",installation_not_implemented, []],
        ["mummer",installation_not_implemented, []],
        ["hmmer",installation_not_implemented, []],
        ["infernal",installation_not_implemented, []],
    )

    for tool in args.skip:
        t.tools[tool].done = True

    t.install(args.tool)

# LOREON installers

Pre-built installers of LOREON **1.1.1** for Windows, macOS and Linux. The files are attached to the [release v1.1.1](https://github.com/iscka/loreon/releases/tag/v1.1.1) on GitHub; the links in the table download them directly.

| System | Installer | Also needed |
|--------|-----------|-------------|
| Windows 10 or 11, 64-bit | [`LOREON_Setup_1.1.1.exe`](https://github.com/iscka/loreon/releases/download/v1.1.1/LOREON_Setup_1.1.1.exe) | Docker Desktop |
| macOS, Apple silicon | [`LOREON-1.1.1-macOS.dmg`](https://github.com/iscka/loreon/releases/download/v1.1.1/LOREON-1.1.1-macOS.dmg) | Python ≥ 3.9, minimap2, samtools |
| Ubuntu 22.04 or 24.04, x86-64 | [`loreon_1.1.1_amd64.deb`](https://github.com/iscka/loreon/releases/download/v1.1.1/loreon_1.1.1_amd64.deb) | minimap2, samtools |

Each installer contains the graphical interface and the pipeline. minimap2 and samtools are not included: on Windows they run inside the Docker image that LOREON builds, on macOS and Linux they are installed as described below.

The installers of the other versions are on the [release page](https://github.com/iscka/loreon/releases). To run LOREON from the source code, see [Installation](../README.md#installation).

## Checking a download

[`SHA256SUMS`](SHA256SUMS) holds the SHA-256 checksum of each installer. Compute the checksum of the downloaded file and compare it with the one in that file:

| System | Command |
|--------|---------|
| Windows (PowerShell) | `Get-FileHash LOREON_Setup_1.1.1.exe` |
| macOS | `shasum -a 256 LOREON-1.1.1-macOS.dmg` |
| Linux | `sha256sum loreon_1.1.1_amd64.deb` |

## Windows

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop) with the WSL 2 backend and start it. Python is not needed.
2. Run `LOREON_Setup_1.1.1.exe`. The installer is not signed, so Windows may warn that its publisher is unknown. If Docker Desktop is not found, the installer says so and offers to open its download page.
3. Start LOREON, tick **Use Docker** and click **Build Image**. The first time, this builds the image in which the analyses run and takes a few minutes; see [Windows / Docker Mode](../README.md#windows--docker-mode).

## macOS

The disk image is for Macs with Apple silicon and macOS 11 or later; it was tested on macOS 27.

1. Install minimap2 and samtools, for example with [Homebrew](https://brew.sh):

   ```bash
   brew install minimap2 samtools
   ```

   Python 3.9 or later is also needed; the one of the Xcode Command Line Tools, which Homebrew needs as well, is suitable.
2. Open `LOREON-1.1.1-macOS.dmg` and drag **LOREON** to the **Applications** folder.
3. Open LOREON. It is not notarized by Apple, so macOS blocks the first start of a downloaded copy. To allow it, open **System Settings > Privacy & Security** and click **Open Anyway**, or run this command in the Terminal and open LOREON again:

   ```bash
   xattr -dr com.apple.quarantine /Applications/LOREON.app
   ```

4. At the first start LOREON creates a Python environment for the pipeline in `~/Library/Application Support/LOREON`. This needs an internet connection and takes one or two minutes.

LOREON looks for minimap2 and samtools in the system folders, in `/opt/homebrew/bin` and in `/usr/local/bin`. If they are somewhere else, for example in a conda environment, start LOREON from a terminal in which they are available:

```bash
open -a LOREON
```

## Linux

The package was built on Ubuntu 22.04 and tested on Ubuntu 22.04 and 24.04.

1. Install the package together with minimap2 and samtools:

   ```bash
   sudo apt install ./loreon_1.1.1_amd64.deb minimap2 samtools
   ```

2. Start LOREON from the applications menu or with the command `loreon`.
3. At the first start LOREON creates a Python environment for the pipeline in `~/.local/share/loreon`. This needs an internet connection and takes one or two minutes.

After the first start, the command `loreon-cli` runs the pipeline from a terminal with the options of the [Command-Line Reference](../README.md#command-line-reference).

The versions of minimap2 and samtools are those of the distribution. LOREON accepts minimap2 2.24 and samtools 1.13 or later; see [System Requirements](../README.md#system-requirements).

## After installing

Click **RUN EXAMPLE DATASET**: LOREON analyses the dataset included in the installer and compares the result with the expected table, as described in [Example Dataset](../README.md#example-dataset).

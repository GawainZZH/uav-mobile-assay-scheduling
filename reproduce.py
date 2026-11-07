"""Copy to the release root as reproduce.py."""
from pathlib import Path
from reproduction40.portable import main

if __name__=='__main__':main(Path(__file__).resolve().parent)

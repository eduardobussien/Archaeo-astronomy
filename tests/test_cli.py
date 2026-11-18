"""
Smoke tests for the command-line tools, run the way a user runs them.

Each tool runs in its own Python process with its output captured through a
pipe, and with that pipe's encoding forced to cp1252, the code page Windows
uses for redirected output. That reproduces on every system the crash the
README's first example once had on Windows ('Δ' has no cp1252 form).
"""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _run(*args):
    env = {k: v for k, v in os.environ.items() if k not in ('PYTHONUTF8', 'PYTHONIOENCODING')}
    env.update({'PYTHONIOENCODING': 'cp1252', 'MPLBACKEND': 'Agg'})
    out = subprocess.run([sys.executable, *args], cwd=ROOT, env=env, capture_output=True)
    assert out.returncode == 0, out.stderr.decode('utf-8', 'replace')
    return out.stdout.decode('utf-8')


def test_the_readme_example_prints_the_star_table_and_the_alignment():
    text = _run('src/alignment.py', '--monument', 'Giza', '--year', '-2559', '--month', '3',
                '--day', '20', '--hour', '21.5')
    assert 'Vondrák 2011' in text
    assert 'Thuban' in text
    assert '(Δ=1.12°)' in text


def test_the_monument_list_names_every_site():
    text = _run('src/alignment.py', '--list-monuments')
    assert text.count('\n  "') == 13
    assert '"Göbekli Tepe"' in text


def test_the_chart_tool_saves_an_image(tmp_path):
    image = tmp_path / 'giza.png'
    text = _run('src/visualize.py', '--monument', 'Giza', '--year', '-2559', '--hour', '21.5',
                '--save', str(image))
    assert 'Saved:' in text
    assert image.read_bytes().startswith(b'\x89PNG')

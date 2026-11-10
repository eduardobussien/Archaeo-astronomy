import os
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _gunicorn_settings(monkeypatch, port):
    if port is None:
        monkeypatch.delenv('PORT', raising=False)
    else:
        monkeypatch.setenv('PORT', port)
    return runpy.run_path(str(ROOT / 'gunicorn.conf.py'))


def test_the_server_listens_on_the_port_the_host_assigns(monkeypatch):
    assert _gunicorn_settings(monkeypatch, '10000')['bind'] == '0.0.0.0:10000'
    assert _gunicorn_settings(monkeypatch, None)['bind'] == '0.0.0.0:8000'


@pytest.mark.parametrize('flask_debug, expected', [(None, False), ('0', False), ('1', True)])
def test_debug_mode_is_off_unless_asked_for(flask_debug, expected):
    env = {k: v for k, v in os.environ.items() if k != 'FLASK_DEBUG'}
    if flask_debug is not None:
        env['FLASK_DEBUG'] = flask_debug
    out = subprocess.run([sys.executable, '-c', 'import app; print(app.app.debug)'],
                         cwd=ROOT, env=env, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == str(expected)

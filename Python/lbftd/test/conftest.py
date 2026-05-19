import os
import pytest

@pytest.fixture(autouse=True)
def _cwd_to_test_dir(monkeypatch):
    monkeypatch.chdir(os.path.dirname(__file__))

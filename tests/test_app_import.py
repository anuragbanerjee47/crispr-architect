import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_app_imports_cleanly():
    app = importlib.import_module("app")
    assert hasattr(app, "main")
    assert callable(app.main)


def test_app_stores_results_in_streamlit_session_state(monkeypatch):
    app = importlib.import_module("app")

    class DummySidebar:
        def header(self, *args, **kwargs):
            return None

        def selectbox(self, label, options, *args, **kwargs):
            if label == "Input Mode":
                return "Paste FASTA Sequence"
            if label == "Cas Enzyme":
                return "SpCas9"
            return kwargs.get("value") or options[0]

        def slider(self, label, *args, **kwargs):
            if "minimum" in label.lower():
                return 40
            if "maximum" in label.lower():
                return 60
            if "mismatch" in label.lower():
                return 3
            return 10

        def text_input(self, label, *args, **kwargs):
            return kwargs.get("value", "")

    class DummyStreamlit:
        def __init__(self):
            self.session_state = {"results": None}
            self.sidebar = DummySidebar()

        def set_page_config(self, *args, **kwargs):
            return None

        def title(self, *args, **kwargs):
            return None

        def button(self, *args, **kwargs):
            return True

        def columns(self, count):
            return [SimpleNamespace(metric=lambda *a, **k: None) for _ in range(count)]

        def subheader(self, *args, **kwargs):
            return None

        def dataframe(self, *args, **kwargs):
            return None

        def table(self, *args, **kwargs):
            return None

        def plotly_chart(self, *args, **kwargs):
            return None

        def warning(self, *args, **kwargs):
            return None

        def download_button(self, *args, **kwargs):
            return None

        def error(self, *args, **kwargs):
            return None

        def selectbox(self, label, options, *args, **kwargs):
            return options[0]

        def text_area(self, *args, **kwargs):
            return "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

    dummy = DummyStreamlit()
    monkeypatch.setattr(app, "st", dummy)

    app.main()

    assert "results" in dummy.session_state
    assert dummy.session_state["results"] is not None
    assert "sequence" in dummy.session_state["results"]
    assert "rows" in dummy.session_state["results"]
    assert "summary" in dummy.session_state["results"]
    assert "vendor_exports" in dummy.session_state["results"]

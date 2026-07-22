"""
Shared headless stubs for PyQt6 / pyvista.

Installed at collection time (before any test module imports
``nics_placer.dialog``) so that whichever test file is collected first
gets the same rich stub set — ``NicsPlacerDialog`` can then be
instantiated and driven directly in tests without a real Qt runtime.

``sys.modules.setdefault`` is used so a test module that wants its own
minimal stub (see ``tests/test_dialog.py``) can still build local stub
classes without clobbering what's already registered here.
"""
import sys
import types
from unittest.mock import MagicMock


class _QBase:
    """Generic stand-in base class for QDialog / QObject.

    Any attribute/method not explicitly defined auto-mocks itself on
    first access (and is cached on the instance), so arbitrary Qt API
    surface (``setWindowTitle``, ``resize``, ``windowFlags``, ...) can
    be called without raising.
    """

    def __init__(self, *args, **kwargs):
        pass

    def __getattr__(self, name):
        val = MagicMock(name=name)
        object.__setattr__(self, name, val)
        return val

    def showEvent(self, event):
        pass

    def closeEvent(self, event):
        pass


class _FakeTableWidgetItem:
    def __init__(self, text=""):
        self._text = text

    def text(self):
        return self._text

    def setText(self, text):
        self._text = text


class _FakeIndex:
    def __init__(self, row):
        self._row = row

    def row(self):
        return self._row


class _FakeTableWidget(_QBase):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._rows = 0
        self._cols = 0
        self._items = {}
        self._selected_rows = []

    def setColumnCount(self, n):
        self._cols = n

    def columnCount(self):
        return self._cols

    def setHorizontalHeaderLabels(self, labels):
        self._header_labels = list(labels)

    def setRowCount(self, n):
        self._rows = n
        self._items = {k: v for k, v in self._items.items() if k[0] < n}

    def rowCount(self):
        return self._rows

    def setItem(self, row, col, item):
        self._items[(row, col)] = item

    def item(self, row, col):
        return self._items.get((row, col))

    def setSelectionBehavior(self, *a):
        pass

    def setEditTriggers(self, *a):
        pass

    def selectedIndexes(self):
        return [_FakeIndex(r) for r in self._selected_rows]

    # test helper (not part of the real QTableWidget API)
    def _set_selected_rows(self, rows):
        self._selected_rows = list(rows)


class _FakeComboBox(_QBase):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._items = []
        self._current = -1
        self.currentIndexChanged = MagicMock()

    def addItem(self, text, data=None):
        self._items.append((text, data))
        if self._current == -1:
            self._current = 0

    def currentData(self):
        if 0 <= self._current < len(self._items):
            return self._items[self._current][1]
        return None

    def findData(self, data):
        for i, (_t, d) in enumerate(self._items):
            if d == data:
                return i
        return -1

    def setCurrentIndex(self, i):
        self._current = i

    def currentIndex(self):
        return self._current

    def blockSignals(self, b):
        return False


def install_qt_stubs():
    """Install rich PyQt6 / pyvista stubs into sys.modules (idempotent)."""
    if "PyQt6" in sys.modules and getattr(sys.modules["PyQt6"], "_nics_placer_rich_stub", False):
        return

    class _Qt:
        class MouseButton:
            LeftButton = 1

        class WindowType:
            WindowMinMaxButtonsHint = 0
            WindowCloseButtonHint = 0

    _QEvent = MagicMock()
    _QEvent.Type.MouseButtonPress = 2
    _QEvent.Type.MouseButtonRelease = 3

    class _QTimer:
        """Not based on _QBase: needs real (non-auto-mocked) default state."""

        _singleshot_calls = []

        def __init__(self, *a, **kw):
            self._interval = None
            self._active = False
            self._timeout_cb = None

        def setInterval(self, ms):
            self._interval = ms

        def isActive(self):
            return self._active

        def start(self, *a):
            self._active = True

        def stop(self):
            self._active = False

        @property
        def timeout(self):
            m = MagicMock()
            m.connect = lambda fn: setattr(self, "_timeout_cb", fn)
            return m

        @staticmethod
        def singleShot(ms, fn):
            _QTimer._singleshot_calls.append((ms, fn))

    class _QCoreApplication:
        _instance = None

        @staticmethod
        def instance():
            return _QCoreApplication._instance

    qt_core = types.ModuleType("PyQt6.QtCore")
    qt_core.Qt = _Qt
    qt_core.QTimer = _QTimer
    qt_core.QObject = _QBase
    qt_core.QEvent = _QEvent
    qt_core.QCoreApplication = _QCoreApplication

    qt_widgets = types.ModuleType("PyQt6.QtWidgets")
    qt_widgets.QDialog = _QBase
    qt_widgets.QAbstractItemView = MagicMock()
    qt_widgets.QComboBox = _FakeComboBox
    qt_widgets.QHBoxLayout = MagicMock()
    qt_widgets.QHeaderView = MagicMock()
    qt_widgets.QLabel = MagicMock()
    qt_widgets.QPushButton = MagicMock()
    qt_widgets.QTableWidget = _FakeTableWidget
    qt_widgets.QTableWidgetItem = _FakeTableWidgetItem
    qt_widgets.QVBoxLayout = MagicMock()

    pyqt6 = types.ModuleType("PyQt6")
    pyqt6.QtCore = qt_core
    pyqt6.QtWidgets = qt_widgets
    pyqt6._nics_placer_rich_stub = True

    for name, mod in [
        ("PyQt6", pyqt6),
        ("PyQt6.QtCore", qt_core),
        ("PyQt6.QtWidgets", qt_widgets),
    ]:
        sys.modules[name] = mod

    class _FakePolyData:
        def __init__(self, positions):
            self.positions = positions
            self.data = {}

        def __setitem__(self, key, value):
            self.data[key] = value

        def glyph(self, geom=None, scale=None, orient=False):
            return MagicMock(name="glyph_mesh")

    pv_module = types.ModuleType("pyvista")
    pv_module.PolyData = _FakePolyData
    pv_module.Sphere = MagicMock(return_value=MagicMock(name="sphere_geom"))
    sys.modules["pyvista"] = pv_module


install_qt_stubs()

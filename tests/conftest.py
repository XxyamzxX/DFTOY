# tests/conftest.py
# -*- coding: utf-8 -*-
"""
Inserta la raíz del proyecto en sys.path para que 'import scf' funcione
al ejecutar pytest desde cualquier carpeta.
POR QUÉ: evita ModuleNotFoundError cuando el paquete no está instalado.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

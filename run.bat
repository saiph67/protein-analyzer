@echo off
rem Start the Protein Analyzer web app (opens http://127.0.0.1:5000 in your browser).
cd /d "%~dp0"
python -m pip install --quiet -r requirements.txt
python app.py
pause

@echo off
rem Builds a stand-alone Windows app in dist\ProteinAnalyzer\ (no Python needed to run it).
rem PyTorch is left out to keep the download small, so the .exe uses the hand-crafted model;
rem the ESM-2 model needs the Python version (python app.py).
cd /d "%~dp0"
python -m pip install --quiet pyinstaller -r requirements.txt || goto :error
python -m PyInstaller --noconfirm --clean --name ProteinAnalyzer --console ^
  --icon NONE ^
  --add-data "templates;templates" --add-data "static;static" --add-data "benchmarks/ss_benchmark.json;benchmarks" ^
  --exclude-module torch --exclude-module transformers --exclude-module tensorflow ^
  --exclude-module IPython --exclude-module tkinter --exclude-module pytest ^
  --distpath dist --workpath build app.py || goto :error
rem Ship the trained results (dataset, hand-crafted model, report, figures) next to the .exe.
robocopy results dist\ProteinAnalyzer\results /E /XF esm_cache_*.joblib model_esm2.joblib jobs.json *.tmp /NFL /NDL /NJH /NJS
if %ERRORLEVEL% GEQ 8 goto :error
echo.
echo Done. Run dist\ProteinAnalyzer\ProteinAnalyzer.exe (or zip the dist\ProteinAnalyzer folder to share it).
exit /b 0
:error
echo Build failed.
exit /b 1

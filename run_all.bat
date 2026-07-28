@echo off
REM Runs every verification script for the Symmetry paper.
setlocal

where python >nul 2>&1
if errorlevel 1 (
  echo Python was not found on PATH. Install Python 3.9 or later and try again.
  exit /b 1
)

echo ============================================================
echo  1/2  verify_all.py
echo ============================================================
python verify_all.py
if errorlevel 1 goto failed

echo.
echo ============================================================
echo  2/2  verify_dedup_independent.py
echo ============================================================
python verify_dedup_independent.py
if errorlevel 1 goto failed

echo.
echo All verification scripts completed successfully.
pause
exit /b 0

:failed
echo.
echo A verification script reported a failure.
pause
exit /b 1

# Check if numpy has been built as a wheel. If it hasn't create one, if it
# has move on.

# NOTE: If you want to force a rebuild of the numpy wheels, modify this file in
# any way and start appveyor again. The cache section in appveyor.yml has set
# this file up as a build dependency.

function InstallNumpy(){

    if (-not(Test-Path "c:\tmp\*.whl")) {
        Write-Host "numpy has not been compiled yet. Starting Long process..."
        Write-Host "pip wheel --wheel-dir=c:\tmp\ numpy"
        iex "cmd /E:ON /V:ON /C .\\appveyor\\run_with_env.cmd pip wheel --wheel-dir=c:\\tmp numpy"
    } else {
        Write-Host "numpy has already been compiled."
        Get-ChildItem "C:\tmp"
    }
}

InstallNumpy

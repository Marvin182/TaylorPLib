- adolc kompilieren:

ben�tigte Komponenten:
Adolc-1.10.2 -> https://projects.coin-or.org/ADOL-C/browser/old/adolc-1.10.2.tar.gz?rev=102
MinGw

MinGw installieren
MinGw Shell starten

$cd /
$mkdir adolc
$mount d:/temp/adolc-1.10.2 /adolc
$cd /adolc/adolc-1.10.2
$./configure
$make
$make install

anschlie�end hat man unter C:\MinGW\msys\1.0\home\<USER> 
den Ordner Adolc-Base

- daeIndexDat kompilieren:

ben�tigte Komponenten:
Dev-C++
MinGw 
TaylorPLibrary (Projekt)

Dev-C++ installieren

den Ordner adolc-base/include/adolc ins Projektverzeichnis kopieren

mit einem Cmd ins Projektvderzeichnis wechseln:

> set PATH=C:\Dev-Cpp\bin\;%PATH%
> C:\Dev-C++\bin\make 
aufrufen

anschlie�end das C:\MinGw\bin zum PATH hinzuf�gen
(;C:\MinGW\bin)

daeIndexDat.exe sollte nun startbar sein

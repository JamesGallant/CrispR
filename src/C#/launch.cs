using System;
using System.Diagnostics;
using System.Windows.Forms;
using System.IO;

/*
https://stackoverflow.com/questions/28174386/how-can-a-bat-file-be-converted-to-exe-without-third-party-tools?answertab=active#tab-top
if csc.exe is not on the path, it can be found in c: windows\microsoft\.NET and can be found there
csc.exe /target:winexe launch.cs /win32icon:icon.ico
The launch.bat must be on the same file location as the launch.cs
*/
class BatCaller {
    static void Main() {
        var batFile = System.Reflection.Assembly.GetEntryAssembly().Location.Replace(".exe", ".cmd");
        if (!File.Exists(batFile)) {
            MessageBox.Show("The launch script could not be found.", "Critical error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            System.Environment.Exit(42);
        }
        var processInfo = new ProcessStartInfo("cmd.exe", "/c \"" + batFile + "\"");
        processInfo.CreateNoWindow = true;
        processInfo.UseShellExecute = false;
        processInfo.RedirectStandardError = true;
        processInfo.RedirectStandardOutput = true;

        var process = Process.Start(processInfo);

        process.OutputDataReceived += (object sender, DataReceivedEventArgs e) => Console.WriteLine("output>>" + e.Data);
        process.BeginOutputReadLine();

        process.ErrorDataReceived += (object sender, DataReceivedEventArgs e) => Console.WriteLine("error>>" + e.Data);
        process.BeginErrorReadLine();

        process.WaitForExit();

        Console.WriteLine("ExitCode: {0}", process.ExitCode);
        process.Close();
    }
}
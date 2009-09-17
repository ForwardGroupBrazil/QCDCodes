//______________________________________________________________________________
void dir(char *path=0)
{
   char s[256] = (!strcmp(gSystem->GetName(), "WinNT")) ? "dir " : "ls -l 
";
   if (path) strcat(s,path);
   gSystem->Exec(s);
}

//______________________________________________________________________________
void lc(char *path=0)
{
   char s[256] = (!strcmp(gSystem->GetName(), "WinNT")) ? "lc " : "ls 
-xstF ";
   if (path) strcat(s,path);
   gSystem->Exec(s);
}


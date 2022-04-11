lfile = "out\x64_Release\Wrapper\Wrapper.lib";
hfile = "include\BS_Simulink.h";
lname = "Wrapper";

clibgen.generateLibraryDefinition ...
    (hfile, "Libraries", lfile, "PackageName", lname);



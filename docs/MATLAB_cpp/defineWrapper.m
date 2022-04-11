%% About defineWrapper.mlx
% This file defines the MATLAB interface to the library |Wrapper|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for &lt;SHAPE&gt;, &lt;DIRECTION&gt;, etc. For more
% information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') Define MATLAB Interface for C++ Library>.



%% Setup. Do not edit this section.
function libDef = defineWrapper()
libDef = clibgen.LibraryDefinition("WrapperData.xml");
%% OutputFolder and Libraries 
libDef.OutputFolder = "C:\Users\kojik\Documents\GitHub\D-BRAIN\docs\MATLAB_cpp";
libDef.Libraries = "out\x64_Release\Wrapper\Wrapper.lib";

%% C++ class |BS_Simulink| with MATLAB name |clib.Wrapper.Simulink| 
SimulinkDefinition = addClass(libDef, "BS_Simulink", "MATLABName", "clib.Wrapper.Simulink", ...
    "Description", "clib.Wrapper.Simulink    Representation of C++ class BS_Simulink."); % Modify help description values as needed.

%% C++ class constructor for C++ class |BS_Simulink| 
% C++ Signature: BS_Simulink::BS_Simulink()
SimulinkConstructor1Definition = addConstructor(SimulinkDefinition, ...
    "BS_Simulink::BS_Simulink()", ...
    "Description", "clib.Wrapper.Simulink    Constructor of C++ class BS_Simulink."); % Modify help description values as needed.
validate(SimulinkConstructor1Definition);

%% C++ class constructor for C++ class |BS_Simulink| 
% C++ Signature: BS_Simulink::BS_Simulink(BS_Simulink const & input1)
SimulinkConstructor2Definition = addConstructor(SimulinkDefinition, ...
    "BS_Simulink::BS_Simulink(BS_Simulink const & input1)", ...
    "Description", "clib.Wrapper.Simulink    Constructor of C++ class BS_Simulink."); % Modify help description values as needed.
defineArgument(SimulinkConstructor2Definition, "input1", "clib.Wrapper.Simulink", "input");
validate(SimulinkConstructor2Definition);

%% C++ function |newSimulink| with MATLAB name |clib.Wrapper.newSimulink|
% C++ Signature: BS_Simulink * newSimulink()
%newSimulinkDefinition = addFunction(libDef, ...
%    "BS_Simulink * newSimulink()", ...
%    "MATLABName", "clib.Wrapper.newSimulink", ...
%    "Description", "clib.Wrapper.newSimulink    Representation of C++ function newSimulink."); % Modify help description values as needed.
%defineOutput(newSimulinkDefinition, "RetVal", "clib.Wrapper.Simulink", <SHAPE>);
%validate(newSimulinkDefinition);

%% C++ function |BS_initialize| with MATLAB name |clib.Wrapper.BS_initialize|
% C++ Signature: void BS_initialize(BS_Simulink * Simulink,int n)
%BS_initializeDefinition = addFunction(libDef, ...
%    "void BS_initialize(BS_Simulink * Simulink,int n)", ...
%    "MATLABName", "clib.Wrapper.BS_initialize", ...
%    "Description", "clib.Wrapper.BS_initialize    Representation of C++ function BS_initialize."); % Modify help description values as needed.
%defineArgument(BS_initializeDefinition, "Simulink", "clib.Wrapper.Simulink", "input", <SHAPE>); % '<MLTYPE>' can be clib.Wrapper.Simulink, or clib.array.Wrapper.Simulink
%defineArgument(BS_initializeDefinition, "n", "int32");
%validate(BS_initializeDefinition);

%% C++ function |BS_sum| with MATLAB name |clib.Wrapper.BS_sum|
% C++ Signature: int BS_sum(BS_Simulink * Simulink)
%BS_sumDefinition = addFunction(libDef, ...
%    "int BS_sum(BS_Simulink * Simulink)", ...
%    "MATLABName", "clib.Wrapper.BS_sum", ...
%    "Description", "clib.Wrapper.BS_sum    Representation of C++ function BS_sum."); % Modify help description values as needed.
%defineArgument(BS_sumDefinition, "Simulink", "clib.Wrapper.Simulink", "input", <SHAPE>); % '<MLTYPE>' can be clib.Wrapper.Simulink, or clib.array.Wrapper.Simulink
%defineOutput(BS_sumDefinition, "RetVal", "int32");
%validate(BS_sumDefinition);

%% Validate the library definition
validate(libDef);

end

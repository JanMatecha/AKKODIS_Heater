model ControledHeaterAuto
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
    Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput SetHeater annotation(
    Placement(transformation(origin = {-120, 10}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor ThermalMass(C = 1000*10, T(start = 293, fixed = true, displayUnit = "K"))  annotation(
    Placement(transformation(origin = {-34, 20}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 10)  annotation(
    Placement(transformation(origin = {-2, 10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput SetTemp annotation(
    Placement(transformation(origin = {80, 12}, extent = {{-20, -20}, {20, 20}}, rotation = 180), iconTransformation(origin = {78, 10}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature prescribedTemperature1 annotation(
    Placement(transformation(origin = {34, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
equation
  connect(SetHeater, prescribedHeatFlow.Q_flow) annotation(
    Line(points = {{-120, 10}, {-80, 10}}, color = {0, 0, 127}));
  connect(prescribedHeatFlow.port, ThermalMass.port) annotation(
    Line(points = {{-60, 10}, {-34, 10}}, color = {191, 0, 0}));
  connect(ThermalMass.port, thermalConductor.port_a) annotation(
    Line(points = {{-34, 10}, {-12, 10}}, color = {191, 0, 0}));
  connect(SetTemp, prescribedTemperature1.T) annotation(
    Line(points = {{80, 12}, {46, 12}}, color = {0, 0, 127}));
  connect(thermalConductor.port_b, prescribedTemperature1.port) annotation(
    Line(points = {{8, 10}, {24, 10}, {24, 12}}, color = {191, 0, 0}));
  annotation(
    uses(Modelica(version = "4.0.0")));
end ControledHeaterAuto;

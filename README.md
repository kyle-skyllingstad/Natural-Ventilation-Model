# Natural-Ventilation-Model

# Overview
This theoretical model predicts temperatures inside Stanford's civil engineering building using a combination of convective and conductive heat flow analysis. It accounts for many factors including internal heat from lights and occupants, night flush convection, radiation, and human-driven airflow. The building, Y2E2, maintains temperatures purely through natural ventilation and does not require any man-made heating or cooling to do so. This bundle predicts these naturally-driven temperatures and compares them to experimental data obtained from sensors mounted to the interior and exterior of the Y2E2 building on Stanford's campus.

# How to Use
In this bundle, there are six Python programs. They are: main host function newBoxModel.py, experimentalData.py, initCondition.py, secondLoop.py, internalHeatFlux.py, and conductionMatrix.py. To run this code and compare theoretical air temperatures to experimental air temperatures, run the main host function newBoxModel.py. No input is needed; this function will call all others and call datasheets from within this chain. It then outputs a visual plot comparing theoretical Y2E2 air temperatures to experimental Y2E2 air temperatures. For more information on what each program doees, see the description at the top of the program.

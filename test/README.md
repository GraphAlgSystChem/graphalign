# Test

This folder contains code to run the program from **src/**. The code in **TestSuite** provides different ways to call the program, either by taking in a string that specifies a specific guide tree or letting the program itself determine the guide tree. Also possible to call the program with and without timeout. 

**main** provides the interface calling the functions provided in **TestSuite**. The scoring function can be redefined by the user. If redefined, ```make``` must be invoked anew.

**GraphIO** is used to load in .json files to create the graphs. See the file for the required format of the json files. 
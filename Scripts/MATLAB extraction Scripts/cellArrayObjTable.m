function [ cellArrayObj ] = cellArrayObjTable(tools)
%Returns table with mapping between (Array Num, Obj Num) coordinates to
%Cell ID


cellArrayObj = [];
cellID = 1; 

while(tools.iterator.continueIteration)
    cellArrayObj = [cellArrayObj; cellID, tools.navigator.currentArrayNum, tools.navigator.currentObjNum];

    tools.iterator.goToNextObject;
    cellID = cellID + 1;


end


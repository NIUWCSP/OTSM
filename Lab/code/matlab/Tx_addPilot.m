function Tx_tilda_addPilot=Tx_addPilot(Tx_tilda,PilotBits);
size_of_Tx_tilda=size(Tx_tilda,1);
Tx_tilda_addPilot=Tx_tilda;
Tx_tilda_addPilot(size_of_Tx_tilda,:)=PilotBits;
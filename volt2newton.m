function volts = volt2newton( volts )
%PrepProd2 convert volt measures into Newtons

volts = volts*0.177918612; %conversion V to KG provided by IT (David) and Katja's linear regression

% volts = volts*0.1667; %conversion V to KG **OLD, before transducer change** 

volts = volts*9.807; %conversion KG to N



end
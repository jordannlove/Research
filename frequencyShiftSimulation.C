/*
 * @author Thomas Gartman
 * @date   September 26th, 2017
 * @file   frequencyShiftSimulation.C
 * @brief  This file should simulate a shift of a frequency from a transmitting point, starting shower point, ending shower point, an ending point
 */

#include <cmath>

void frequencyShiftSimulation()
{
	//variables
	//Distance Units are in Meters
	//Angles are in degrees
	//time is in seconds
	Double_t transmitterLocationX = 0;
	Double_t transmitterLocationY = 0;
	Double_t transmitterLocationZ = 0;
	Double_t receiverLocationX = 0;
	Double_t receiverLocationY = 0;
	Double_t receiverLocationZ = 0;
	Double_t theta = 0;
	Double_t phi = 0;
        Double_t showerStartingX = 0;
	Double_t showerStartingY = 0;
	Double_t showerStartingZ = 0;
	Double_t showerEndingX = 0;
	Double_t showerEndingY = 0;
	Double_t showerEndingZ = 0;
	Double_t length0 = 0;
	Double_t length1 = 0;
	Double_t length2 = 0;
	Double_t length3 = 0;
	Double_t beta1 = 0;
	Double_t beta2 = 0;
	Double_t frequencyShift1 = 0;
	Double_t frequencyShift2 = 0;
	
	const Double_t DISTANCESHOWER = 30.0;
	const Double_t SPEEDOFLIGHT = 300000000.0;
	const Double_t TIMESHOWER = 30.0/300000000.0;
 
	TH1D* grShifts = new TH1D("FinalFrequency Graph", "FinalFrequency Vs Event Number", 100000, 1, 100000);
        TH1D* histShifts = new TH1D("FinalFrequency Histogram", "FinalFrequency Histogram", 61, -50, 3000);	
	TDatime* time = new TDatime();

	Double_t RxArray[100000];
	Double_t TxArray[100000];
	Double_t TxRxArray[100000];
        Double_t thetaArray[100000];
	Double_t phiArray[100000];
	Double_t frequencyShiftArray[100000];
	TRandom2* eventGenerator = new TRandom2(time->Get());
	
	TFile* outFile = new TFile("FrequencyShiftSimulation.root", "recreate");

	//randomize Locations
        transmitterLocationX = eventGenerator->Uniform(-2000,2000);
	transmitterLocationY = eventGenerator->Uniform(-2000,2000);
	transmitterLocationZ = eventGenerator->Uniform(-2000,2000);
	receiverLocationX = eventGenerator->Uniform(-2000,2000);
	receiverLocationY = eventGenerator->Uniform(-2000,2000);
	receiverLocationZ = eventGenerator->Uniform(-2000,2000);
	for(int i = 1; i < 100000; i++)
	{
	  theta = eventGenerator->Uniform(0, 360);
	  phi = eventGenerator->Uniform(0, 360);
	  showerStartingX = eventGenerator->Uniform(-2000,2000);
	  showerStartingY = eventGenerator->Uniform(-2000,2000);
	  showerStartingZ = eventGenerator->Uniform(-2000,2000);

	  //calculate ending shower location
	  showerEndingX = DISTANCESHOWER * sin(theta) * cos(phi) + showerStartingX;
	  showerEndingY = DISTANCESHOWER * sin(theta) * sin(phi) + showerStartingY;
	  showerEndingZ = DISTANCESHOWER * cos(theta) + showerStartingZ;

	  //Calculate dL and the frequency shift
	  //Calculate the shift wrt to the receiver and the shift wrt to the transmitter
          length0 = sqrt(pow(showerStartingX - transmitterLocationX, 2) + pow(showerStartingY - transmitterLocationY, 2) + pow(showerStartingZ - transmitterLocationZ, 2));
	  length1 = sqrt(pow(showerEndingX - transmitterLocationX, 2) + pow(showerEndingY - transmitterLocationY, 2) + pow(showerEndingZ - transmitterLocationZ, 2));
	  length2 = sqrt(pow(showerStartingX - receiverLocationX, 2) + pow(showerStartingY - receiverLocationY, 2) + pow(showerStartingZ - receiverLocationZ, 2)); 
	  length3 = sqrt(pow(showerEndingX - receiverLocationX, 2) + pow(showerEndingY - receiverLocationY, 2) + pow(showerEndingZ - receiverLocationZ, 2));

	  beta1 = ((length1 - length0)/TIMESHOWER)/SPEEDOFLIGHT;
	  beta2 = ((length3 - length2)/TIMESHOWER)/SPEEDOFLIGHT;
	  
	  TxArray[i-1] = beta1;
	  RxArray[i-1] = beta2;
	  TxRxArray[i-1] = beta1*beta2;
          thetaArray[i-1] = theta;
	  phiArray[i-1] = phi;

	  frequencyShift1 = 450 * sqrt((1-beta1)/(1+beta1));
          frequencyShift2 = frequencyShift1 * sqrt((1+beta2)/(1-beta2));
	  
	  frequencyShiftArray[i] = frequencyShift2;

	  grShifts->AddBinContent(i, frequencyShift2);
	}
	TH2D* gr2DDoppler = new TH2D("gr2DDoppler", "TxVsRxGraph", 200, -1, 1, 200, -1, 1);
	gr2DDoppler->SetTitle("Tx vs Rx Graph");
	TH2D* grRxTheta = new TH2D("TxRx", "TxRxTheta", 360, 0, 360, 200, -1, 1);
	grRxTheta->SetTitle("TxRx Vs Theta");
	TH2D* grRxPhi = new TH2D("TxRxPhi", "TxRxPhi", 360, 0, 360, 200, -1, 1);
	grRxPhi->SetTitle("TxRx Vs Phi");

        histShifts->SetTitle("Final Frequency Histogram");

	for(int i = 0; i <= 100000; i++)
	{
	  gr2DDoppler->Fill(TxArray[i], RxArray[i]);
	  grRxTheta->Fill(thetaArray[i], TxRxArray[i]);
	  grRxPhi->Fill(phiArray[i], TxRxArray[i]);
	  histShifts->Fill(frequencyShiftArray[i], 1);
	  
	}

	grRxPhi->Write("TxRxVsPhi");
	gr2DDoppler->Write("TxVsRxGraph");
	grRxTheta->Write("TxRxVsTheta");
	grShifts->Write("FinalFrequencyVsEventNumber");
        histShifts->Write("FinalFrequencyHistogram");

	outFile->Close();
	delete grShifts;
	delete time;
	delete eventGenerator;
	delete outFile;
}

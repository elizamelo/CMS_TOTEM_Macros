import itertools
from ROOT import TCanvas,TFile, TH1, TSlider,TLegend,TString,TLatex
from ROOT import gROOT, gBenchmark, gRandom,gStyle

#https://svnweb.cern.ch/trac/penn/browser/reece/rel/tree_trimmer/trunk/tree_trimmer.py
gROOT.Reset();
gROOT.ProcessLine('.L styleTDR.C')
gROOT.ProcessLine('setTDRStyle()')
#gROOT.ForceStyle() 
#gStyle.SetOptStat(0)
gBenchmark.Start( 'Plotter' )
# Create a new canvas, and customize it.
#c1 = TCanvas( 'c1', 'Control Plots', 200, 10, 600, 400 )
#c1.SetGrid();


## List of file
filesData = 'Data_Mll_6_10_RPPlusAccept.root'
#filesData = 'Data_Mll_10_20_RPPlusAccept.root'
#filesPompty = ['POMPTY_Mll_6_10_RPPlusAccept.root','POMPTY_Mll_10_20_RPPlusAccept.root','POMPTY_Mll_20_60_RPPlusAccept.root'] 
#filesPythia = ['PYTHIA_Mll_6_10_RPPlusAccept.root','PYTHIA_Mll_10_20_RPPlusAccept.root','PYTHIA_Mll_20_60_RPPlusAccept.root'] 
filesPompty = ['POMPTY_Mll_6_10_RPPlusAccept.root']#,'POMPTY_Mll_10_20_RPPlusAccept.root','POMPTY_Mll_20_60_RPPlusAccept.root']
filesPythia = ['PYTHIA_Mll_6_10_RPPlusAccept.root']#,'PYTHIA_Mll_10_20_RPPlusAccept.root','PYTHIA_Mll_20_60_RPPlusAccept.root']
#filesData = 

#data_hist_name = 

hist_name = ['dy_mass','dy_eta','dy_pt','muon1_pt','muon1_eta' ]
X_axis=["M_{#mu^{+}#mu^{-}} [GeV]","#eta_{#mu^{+}#mu^{-}}","p_{T_{#mu^{+}#mu^{-}}} [GeV]","p_{T_{#mu^{1}}} [GeV]","#eta_{#mu^{1}}"]


#hist_name = ['dy_mass']
#data_hist_name = 'dimuon_mass_proton_left'

#hist_name = 'dy_mass_rpplus_accept'
#data
#XSEC M_6_10 Plus
#BR= 0.04 #4%
BR =1.0
xsec_pompty = 7.70E-1 # units in nb
xsec_pythia = 7.48 # units in nb
#Number of Events M_6_10 Plus
n_gen_pompty = 708224
n_gen_pythia = 92500 ## Total GEN Number 
#n_gen_pythia = 13006 ## in DY region 
LumiEff = 49.15 # nb-1
Trigger_eff = 1.0

weight_data = 1.0/LumiEff

#filesData = [   ]

#pompty
weight_pompty = (xsec_pompty)/(n_gen_pompty)
#weight_pompty = ((xsec_pompty*BR) * LumiEff)/(n_gen_pompty*Trigger_eff)
# pythia
weight_pythia = (xsec_pythia)/(n_gen_pythia)
#weight_pythia = ((xsec_pythia*BR) * LumiEff)/(n_gen_pythia*Trigger_eff)

#print ' weight_pompty: ', weight_pompty
#print ' weight_pythia: ', weight_pythia

filesData=TFile.Open( filesData, 'READ' )

#####################
## Loop for histo
for h_name,xaxis in itertools.izip (hist_name,X_axis):
#for h_name in hist_name:
    #for xaxis in X_axis:      		
	#hfile = TFile( h_name + '.root', 'RECREATE', 'Demo ROOT file with histograms' )
	# Canvas
	c1 = TCanvas( 'c1', 'Control Plots', 200, 10, 600, 400 )
  	gStyle.SetOptStat(0)
	gROOT.ForceStyle()
	#Legend
 	leg = TLegend(0.6946309,0.7446237,0.9395973,0.922043);
 	leg.SetBorderSize(0);
 	leg.SetLineStyle(0);
 	leg.SetTextFont(42);
 	leg.SetFillStyle(0);
 	leg.SetFillColor(0);

	# INPUT DATA
#	filesData=TFile.Open( filesData, 'READ' )
	#print ('Data Histograms: ',h_name,xaxis)
	h_data = filesData.Get(h_name)
#	h_data.Sumw2()
	h_data.SetStats(False)
	h_data.SetMarkerStyle(20)
	h_data.SetMarkerSize(1.3)
	h_data.GetYaxis().SetTitle("N Events")
  	h_data.GetXaxis().SetTitle(xaxis)
	#TString titlex = h_data.GetXaxis().GetTitle();
	leg.AddEntry(h_data,"Data 2012 Low PU ","P")
	
	##Loop for MC files
        for pompty_name,pythia_name in itertools.izip(filesPompty,filesPythia):
#	for pompty_name in filesPompty:
		file_pompty=TFile.Open( pompty_name, 'READ' )
		#file_pompty.ls()
		print ('POMPTY Histograms: ',h_name,xaxis) 
		h_pompty = file_pompty.Get(h_name)
#		h_pompty.Sumw2()
	#    	h_pompty.Scale(1.0/h_pompty.GetEntries())
		h_pompty.SetStats(False)
		h_pompty.SetFillColor(2);
      		h_pompty.SetLineWidth(2);
        	h_pompty.SetFillStyle(3004);
        	h_pompty.GetYaxis().SetTitle("N Events")
		h_pompty.GetXaxis().SetTitle(xaxis)
		leg.AddEntry(h_pompty,"Pompty Plus","LFP")
#		h_pompty.Draw('histosame')

#	for pythia_name in filesPythia:
		file_pythia=TFile.Open( pythia_name, 'READ' )
		# file_pythia.ls()
		print ('PYTHIA Histograms: ',h_name,xaxis)	
		h_pythia = file_pythia.Get(h_name)
#		h_pythia.Sumw2()
		#h_pythia.Scale(1.0/h_pythia.GetEntries())
 		h_pythia.SetStats(False)	
		h_pythia.SetFillColor(4)
        	h_pythia.SetLineWidth(2)
 		h_pythia.SetFillStyle(3004)
 		h_pythia.GetYaxis().SetTitle("N Events")
		h_pythia.GetXaxis().SetTitle(xaxis)
		leg.AddEntry(h_pythia,"PYTHIA DY","LFP")
	#	h_pythia.Draw('histosame')
	#	h_pythia.Draw('histo')
############################################
#Normalization per events
#############################################
	
	weight0 = 1./h_data.GetEntries();#DATA
	weight1 = 1./h_pompty.GetEntries();#POMPTY
	weight2 = 1./h_pythia.GetEntries();# PYTHIA
 

	ratio_pompty = weight1/weight0 # MC/DATA (MC dist. UP)
	h_pompty.Scale(ratio_pompty);

	ratio_pythia = weight2/weight0 # MC/DATA (MC dist. UP)

	h_pythia.Scale(ratio_pythia);
	print "normalization for n of evts"
##############
## http://stackoverflow.com/questions/1748641/boolean-in-python

##		h_pompty.Scale(weight_pompty)
##		h_pythia.Scale(weight_pythia)
#################################################
	h_data.Draw('')
	h_pompty.Draw('histosame')
	h_pythia.Draw('histosame')
	leg.Draw("same") 
	c1.Update();
	# Text File
	text1 = TLatex(3.570061,23.08044,"Work in Progress")
	text1.SetNDC()
	text1.SetTextAlign(13)
	text1.SetX(0.184)
	text1.SetY(0.928)
#         //text1->SetLineWidth(2);
	text1.SetTextFont(42)
	text1.SetTextSizePixels(24)
	text1.Draw()
	#hfile.Write()
	c1.Modified()
	c1.Update()
	c1.SaveAs(h_name + '.root')
#	c1.Modified()
	##Save histogram
	#hfile = TFile( h_name + '.root', 'RECREATE', 'Demo ROOT file with histograms' )
	#h1f.Write()
'''
# Create some histograms.
total  = TH1F( 'total', 'This is the total distribution', 100, -4, 4 )
main   = TH1F( 'main', 'Main contributor', 100, -4, 4 )
s1     = TH1F( 's1', 'This is the first signal', 100, -4, 4 )
s2     = TH1F( 's2', 'This is the second signal', 100, -4, 4 )
total.Sumw2()   # this makes sure that the sum of squares of weights will be stored

# Set canvas/frame attributes.
	total.SetMarkerStyle( 21 )
	total.SetMarkerSize( 0.7 )
	main.SetFillColor( 16 )
	s1.SetFillColor( 42 )
s2.SetFillColor( 46 )

# Initialize random number generator.
gRandom.SetSeed()
	gauss, landau = gRandom.Gaus, gRandom.Landau

# for speed, bind and cache the Fill member functions
	histos = [ 'total', 'main', 's1', 's2' ]
	for name in histos:
	exec '%sFill = %s.Fill' % (name,name)

# Fill histograms randomly
	kUPDATE = 500
	for i in xrange( 10000 ):
# Generate random values.
	xmain = gauss( -1, 1.5 )
	xs1   = gauss( -0.5, 0.5 )
	xs2   = landau( 1, 0.15 )
mainFill( xmain )

# Fill histograms.
	s1Fill( xs1, 0.3 )
	s2Fill( xs2, 0.2 )
	totalFill( xmain )
	totalFill( xs1, 0.3 )
totalFill( xs2, 0.2 )

# Update display every kUPDATE events.
	if i and (i%kUPDATE) == 0 :
	if i == kUPDATE :
	total.Draw( 'e1p' )
	main.Draw( 'same' )
	s1.Draw( 'same' )
	s2.Draw( 'same' )
c1.Update()
	slider = TSlider( 'slider', 'test', 4.2, 0, 4.6, total.GetMaximum(), 38 )
slider.SetFillColor( 46 )

	if slider:
slider.SetRange( 0, float(i) / 10000. )

	c1.Modified()
c1.Update()

# Destroy member functions cache.
	for name in histos:
	exec 'del %sFill' % name
	del histos

# Done, finalized and trigger an update.
slider.SetRange( 0, 1 )
	total.Draw( 'sameaxis' ) # to redraw axis hidden by the fill area
	c1.Modified()
c1.Update()
	'''
gBenchmark.Show( 'Plotter' )


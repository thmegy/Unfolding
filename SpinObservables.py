import ROOT
import math

file = ROOT.TFile('/home/theo/Téléchargements/ttbar_dilep_0j_SM.root', 'READ')
tree = file.Get('truth')

h = ROOT.TH1F('', '', 20, -1, 1)


# boucle sur les événements d'un TTree
for i, evt in enumerate(tree):

    if (i%10000==0):
        print(i)
    #Find ttbar rest frame in order to define the spin quantisation axes

    # Define top and anti-top in lab frame

    t = ROOT.TLorentzVector()
    tbar = ROOT.TLorentzVector()
    t.SetPtEtaPhiM(evt.t_pt, evt.t_eta, evt.t_phi, evt.t_m)
    tbar.SetPtEtaPhiM(evt.tbar_pt, evt.tbar_eta, evt.tbar_phi, evt.tbar_m)

    # Build ttbar system and find boost needed to go to its rest frame (spatial components divided by time component)
    ttbar = t + tbar
    boost_to_ttbar_CM = ttbar.BoostVector()

    # Move top to ttbar rest frame ---> Get spin quantisation axis k
    t.Boost(-boost_to_ttbar_CM)
    tbar.Boost(-boost_to_ttbar_CM)

    #Find boosts to respective top and anti-top rest frames (spatial components divided by time component)
    boost_to_t_CM = t.BoostVector()
    boost_to_tbar_CM = tbar.BoostVector()



    k_axis = t.Vect()
    k_axis *= 1/k_axis.Mag()  # normalise vector to 1

    # Define beam axis direction ---> Get spin quantisation axes n and r
    beam_axis = ROOT.TVector3()
    beam_axis.SetXYZ(0,0,1)

    y = k_axis.Dot(beam_axis)
    r = math.sqrt(1-y*y)

    n_axis = beam_axis.Cross(k_axis)
    n_axis *= 1/r
    n_axis *= y/abs(y)  # multiply by sign of y                                                       

    r_axis = (beam_axis-y*k_axis)
    r_axis *= 1/r
    r_axis *= y/abs(y)  # multiply by sign of y 

    
    # Get angles needed for spin polarisations and correlations

    # Define lplus and lminus in lab frame, then move them to top and anti-top rest frames, respectively
    lplus = ROOT.TLorentzVector()
    lminus = ROOT.TLorentzVector()

    lplus.SetPtEtaPhiM(evt.lbar_pt, evt.lbar_eta, evt.lbar_phi, evt.lbar_m)
    lminus.SetPtEtaPhiM(evt.l_pt, evt.l_eta, evt.l_phi, evt.l_m)
      
    cos_theta_k_plus = 0

    lplus.Boost(-boost_to_ttbar_CM)
    lminus.Boost(-boost_to_ttbar_CM)
    lplus.Boost(-boost_to_t_CM)
    lminus.Boost(-boost_to_tbar_CM)

    cos_theta_k_plus = math.cos(lplus.Vect().Angle(k_axis))
    cos_theta_k_minus = math.cos(lminus.Vect().Angle(-k_axis))

    corr_kk = cos_theta_k_plus*cos_theta_k_minus

    h.Fill(corr_kk)


h.Draw()
print(-9*h.GetMean())
input('')

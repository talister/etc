import os
from astropy.config import ConfigNamespace, ConfigItem

__all__ = ['conf']


class Conf(ConfigNamespace):
    """Configuration parameters."""

    bessell_U_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.U', 'Bessell U')
    bessell_B_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.B', 'Bessell B')
    bessell_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.V', 'Bessell V')
    bessell_R_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.R', 'Bessell R')
    bessell_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Generic/Bessell.I', 'Bessell I')

    lco_u_file = ConfigItem('comp/lco/SDSS.up.txt', 'LCO SDSS u')
    lco_g_file = ConfigItem('comp/lco/SDSS.gp.txt', 'LCO SDSS g')
    lco_r_file = ConfigItem('comp/lco/SDSS.rp.txt', 'LCO SDSS r')
    lco_i_file = ConfigItem('comp/lco/SDSS.ip.txt', 'LCO SDSS i')
    lco_zs_file = ConfigItem('comp/lco/PSTR-ZS-avg.txt', 'LCO SDSS/PanSTARRS zs')
    lco_w_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=LasCumbres/LasCumbres.PS_w', 'LCO PanSTARRS w')
    lco_Y_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=LasCumbres/LasCumbres.PS_y', 'LCO PanSTARRS Y')

    lco_c2_file = ConfigItem("comp/lco/LCO_ESA_C2.csv", "LCO ESA C2")
    lco_c3_file = ConfigItem("comp/lco/LCO_ESA_C3.csv", "LCO ESA C3")
    lco_oh_file = ConfigItem("comp/lco/LCO_ESA_OH.csv", "LCO ESA OH")
    lco_cn_file = ConfigItem("comp/lco/LCO_ESA_CN.csv", "LCO ESA CN")
    lco_nh2_file = ConfigItem("comp/lco/LCO_ESA_NH2.csv", "LCO ESA NH2")
    lco_cr_file = ConfigItem("comp/lco/LCO_ESA_CR.csv", "LCO ESA CR")

    lco_fli_clear_file = ConfigItem("comp/lco/LCO_ADNP-P1-002.csv", "LCO Astrodon clear")
    lco_fli_U_file = ConfigItem("comp/lco/LCO_ADNP-UV-001.csv", "LCO Astrodon U")
    lco_fli_B_file = ConfigItem("comp/lco/LCO_ADNP-BU-001.csv", "LCO Astrodon B")
    lco_fli_V_file = ConfigItem("comp/lco/LCO_ADNP-VX-001.csv", "LCO Astrodon V")
    lco_fli_R_file = ConfigItem("comp/lco/LCO_ADNP-RS-001.csv", "LCO Astrodon R")
    lco_fli_I_file = ConfigItem("comp/lco/LCO_ADNP-IC-001.csv", "LCO Astrodon I")
    lco_fli_g_file = lco_g_file
    lco_fli_r_file = lco_r_file
    lco_fli_i_file = lco_i_file
    lco_fli_zs_file = lco_zs_file

    lco_U_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-ux.txt', 'LCO Bessell U')
    lco_B_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-bx.txt', 'LCO Bessell B')
    lco_V_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-vx.txt', 'LCO Bessell V')
    lco_R_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-rx.txt', 'LCO Bessell R')
    lco_I_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-ix.txt', 'LCO Bessell I')

    wht_U_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.RGO_U9', 'WHT/PFIP RGO U9')
    wht_B_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_B', 'WHT/PFIP Harris B')
    wht_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_V', 'WHT/PFIP Harris V')
    wht_R_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_R', 'WHT/PFIP Harris R')
# This abruptly stops at 900nm and ~82% transmission causing issues.
#    wht_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_I', 'WHT/PFIP Harris I')
# This second version is just the glass from Schott, scaled from 3->4mm but
# reproduces above very well and continues down to 0 transmission and so is
# better behaved.
    wht_I_file = ConfigItem("comp/ing/WHT_Harris_I.dat", 'WHT/PFIP Harris I')

    ctio_U_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.bessel_U', 'CTIO/SOI Bessell U')
    ctio_B_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.bessel_B', 'CTIO/SOI Bessell B')
    ctio_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.bessel_V', 'CTIO/SOI Bessell V')
    ctio_R_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.bessel_R', 'CTIO/SOI Bessell R')
    ctio_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.bessel_I', 'CTIO/SOI Bessell I')

    ctio_u_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.sdss_u', 'CTIO/SOI SDSS u')
    ctio_g_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.sdss_g', 'CTIO/SOI SDSS g')
    ctio_r_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.sdss_r', 'CTIO/SOI SDSS r')
    ctio_i_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.sdss_i', 'CTIO/SOI SDSS i')
    ctio_z_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.sdss_z', 'CTIO/SOI SDSS z')
    ctio_CN_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=CTIO/SOI.CN', 'CTIO/SOI CN')

    soar_goodman_600l_file = ConfigItem('comp/soar/Goodman_600l_grating.dat', 'SOAR Goodman 600lines/mm grating order=1')

    lsst_u_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_u.dat', 'LSST u')
    lsst_g_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_g.dat', 'LSST g')
    lsst_r_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_r.dat', 'LSST r')
    lsst_i_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_i.dat', 'LSST i')
    lsst_z_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_z.dat', 'LSST z')
    lsst_y_file = ConfigItem('$CDBS_PATH/comp/lsst/filter_y.dat', 'LSST y')

    eso_U_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.0877', 'ESO/WFI U')
    eso_B_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.0878', 'ESO/WFI B')
    eso_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.0843', 'ESO/WFI V')
    eso_Rc_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.0844', 'ESO/WFI Rc')
    eso_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.0879', 'ESO/WFI I')

    eso_uHIGH_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1112', 'ESO/FORS u_HIGH')
    eso_bHIGH_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1113', 'ESO/FORS b_HIGH')
    eso_vHIGH_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1114', 'ESO/FORS v_HIGH')
    eso_gHIGH_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1115', 'ESO/FORS g_HIGH')
    eso_fors_R_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1076', 'ESO/FORS Bessell R Special')
    eso_fors_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1077', 'ESO/FORS Bessell I')
    eso_fors_z_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=ESO/ESO.1078', 'ESO/FORS Gunn z')
    eso_vst_uprime_file = ConfigItem('comp/eso/sloan_u_prime.dat', "ESO Omegacam u'")
    eso_vst_gprime_file = ConfigItem('comp/eso/sloan_g_prime.dat', "ESO Omegacam g'")
    eso_vst_rprime_file = ConfigItem('comp/eso/sloan_r_prime.dat', "ESO Omegacam r'")
    eso_vst_iprime_file = ConfigItem('comp/eso/sloan_i_prime.dat', "ESO Omegacam i'")
    eso_vst_zprime_file = ConfigItem('comp/eso/sloan_z_prime.dat', "ESO Omegacam z'")
    # These are the more correct versions, having being extracted from the 47
    # point scans across the filter and averaged but the above, which if I am
    # correct, include the CCD response, is what are used in the ESO ETC and then
    # partly undone through use of the Omegacam optics fudge file... <sigh>
    # eso_vst_u_file = ConfigItem('comp/eso/ESO_Omegacam_u.csv', "ESO Omegacam u'")
    # eso_vst_g_file = ConfigItem('comp/eso/ESO_Omegacam_g.csv', "ESO Omegacam g'")
    # eso_vst_r_file = ConfigItem('comp/eso/ESO_Omegacam_r.csv', "ESO Omegacam r'")
    # eso_vst_i_file = ConfigItem('comp/eso/ESO_Omegacam_i.csv', "ESO Omegacam i'")
    # eso_vst_z_file = ConfigItem('comp/eso/ESO_Omegacam_z.csv', "ESO Omegacam z'")

    eso_fors_300V_file = ConfigItem("comp/eso/ESO_FORS2_Grism_300V.dat")
    eso_fors_600B_file = ConfigItem("comp/eso/ESO_FORS2_Grism_600B.dat")
    eso_efosc2_grism1_file = ConfigItem("comp/eso/ESO_EFOSC2_Grism1.dat")
    eso_efosc2_grism2_file = ConfigItem("comp/eso/ESO_EFOSC2_Grism2.dat")
    eso_feros_file = ConfigItem("comp/eso/ESO_FEROS_slitfiber.dat")

    gemini_gmosn_g_file = ConfigItem("comp/gemini/gmos_n_g_G0301.txt", 'Gemini/GMOS-N g')
    gemini_gmosn_r_file = ConfigItem("comp/gemini/gmos_n_r_G0303.txt", 'Gemini/GMOS-N r')
# These don't go down to 0 at the ends of the bandpass, creating extrapolation problems
#    gemini_gmosn_g_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Gemini/GMOS-N.g', 'Gemini/GMOS-N g')
#    gemini_gmosn_r_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Gemini/GMOS-N.r', 'Gemini/GMOS-N r')
    gemini_gmosn_i_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Gemini/GMOS-N.i', 'Gemini/GMOS-N i')
    gemini_gmosn_z_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Gemini/GMOS-N.z', 'Gemini/GMOS-N z')
    gemini_gmosn_ri_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=Gemini/GMOS-N.ri', 'Gemini/GMOS-N ri')

    optics_NaCl_file = ConfigItem("comp/optics/NaCl.dat", "NaCl")
    optics_BAK2_file = ConfigItem("comp/optics/BAK2_glass.dat", "BAK2 glass")
    optics_CaF2_file = ConfigItem("comp/optics/CaF2_glass.dat", "CaF2 glass")
    optics_FK5_file = ConfigItem("comp/optics/FK5_glass.dat", "FK5 crown glass")
    optics_FK58_file = ConfigItem("comp/optics/FK58_glass.dat", "FK58 crown glass")
    optics_LAL7_file = ConfigItem("comp/optics/LAL7_glass.dat", "LAL7 glass")
    optics_UVFS_file = ConfigItem("comp/optics/UVFusedSilica.dat", "UV Fused Silica")

    mapping = {
                'U' : bessell_U_file,
                'B' : bessell_B_file,
                'V' : bessell_V_file,
                'R' : bessell_R_file,
                'I' : bessell_I_file,
                'u' : lco_u_file,
                'g' : lco_g_file,
                'r' : lco_r_file,
                'i' : lco_i_file,
                'z' : lco_zs_file,
                'zs' : lco_zs_file,
                'w'  : lco_w_file,
                'Y'  : lco_Y_file,
                'C2' : lco_c2_file,
                'C3' : lco_c3_file,
                'OH' : lco_oh_file,
                'CN' : lco_cn_file,
                'NH2': lco_nh2_file,
                'CR' : lco_cr_file,
                'LCO::U' : lco_U_file,
                'LCO::B' : lco_B_file,
                'LCO::V' : lco_V_file,
                'LCO::R' : lco_R_file,
                'LCO::I' : lco_I_file,
                'LCO::FLI::clear' : lco_fli_clear_file,
                'LCO::FLI::U' : lco_fli_U_file,
                'LCO::FLI::B' : lco_fli_B_file,
                'LCO::FLI::V' : lco_fli_V_file,
                'LCO::FLI::R' : lco_fli_R_file,
                'LCO::FLI::I' : lco_fli_I_file,
                'LCO::FLI::gp' : lco_fli_g_file,
                'LCO::FLI::rp' : lco_fli_r_file,
                'LCO::FLI::ip' : lco_fli_i_file,
                'LCO::FLI::zs' : lco_fli_zs_file,
                'WHT::U' : wht_U_file,
                'WHT::B' : wht_B_file,
                'WHT::V' : wht_V_file,
                'WHT::R' : wht_R_file,
                'WHT::I' : wht_I_file,
                'CTIO::U' : ctio_U_file,
                'CTIO::B' : ctio_B_file,
                'CTIO::V' : ctio_V_file,
                'CTIO::R' : ctio_R_file,
                'CTIO::I' : ctio_I_file,
                'CTIO::u' : ctio_u_file,
                'CTIO::g' : ctio_g_file,
                'CTIO::r' : ctio_r_file,
                'CTIO::i' : ctio_i_file,
                'CTIO::z' : ctio_z_file,
                'CTIO::CN' : ctio_CN_file,
                'SOAR:600l/mm' : soar_goodman_600l_file,
                'ESO::U' : eso_U_file,
                'ESO::B' : eso_B_file,
                'ESO::V' : eso_V_file,
                'ESO::Rc' : eso_Rc_file,
                'ESO::I' : eso_I_file,
                'ESO::Omegacam::up' : eso_vst_uprime_file,
                'ESO::Omegacam::gp' : eso_vst_gprime_file,
                'ESO::Omegacam::rp' : eso_vst_rprime_file,
                'ESO::Omegacam::ip' : eso_vst_iprime_file,
                'ESO::Omegacam::zp' : eso_vst_zprime_file,
                # 'ESO::Omegacam::u' : eso_vst_u_file,
                # 'ESO::Omegacam::g' : eso_vst_g_file,
                # 'ESO::Omegacam::r' : eso_vst_r_file,
                # 'ESO::Omegacam::i' : eso_vst_i_file,
                # 'ESO::Omegacam::z' : eso_vst_z_file,
                'ESO::FORS::u' : eso_uHIGH_file,
                'ESO::FORS::b' : eso_bHIGH_file,
                'ESO::FORS::v' : eso_vHIGH_file,
                'ESO::FORS::g' : eso_gHIGH_file,
                'ESO::FORS::R' : eso_fors_R_file,
                'ESO::FORS::I' : eso_fors_I_file,
                'ESO::FORS::z' : eso_fors_z_file,
                'ESO::FORS::300V' : eso_fors_300V_file,
                'ESO::FORS::600B' : eso_fors_600B_file,
                'ESO::EFOSC2::Grism1' : eso_efosc2_grism1_file,
                'ESO::EFOSC2::Grism2' : eso_efosc2_grism2_file,
                'ESO::FEROS' : eso_feros_file,
                'LSST::u' : lsst_u_file,
                'LSST::g' : lsst_g_file,
                'LSST::r' : lsst_r_file,
                'LSST::i' : lsst_i_file,
                'LSST::z' : lsst_z_file,
                'LSST::y' : lsst_y_file,
                'Gemini::GMOS-N::g' : gemini_gmosn_g_file,
                'Gemini::GMOS-N::r' : gemini_gmosn_r_file,
                'Gemini::GMOS-N::i' : gemini_gmosn_i_file,
                'Gemini::GMOS-N::z' : gemini_gmosn_z_file,
                'Gemini::GMOS-N::ri' : gemini_gmosn_ri_file,
                'NaCl' : optics_NaCl_file,
                'UVFS' : optics_UVFS_file,
                'CaF2' : optics_CaF2_file,
                'BAK2' : optics_BAK2_file,
                'FK5'  : optics_FK5_file,
                'FK58'  : optics_FK58_file,
                'LAL7'  : optics_LAL7_file,
              }

    # STANDARD STARS
    vega_file = ConfigItem(
        'http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_010.fits', 'Vega')
    sun_file = ConfigItem(os.path.join('$CDBS_PATH', 'calspec', 'sun_reference_stis_002.fits'), "Solar reference spectrum from https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec.html")

    # These are generated from a modified version of Lynne Jones's astcolors.py
    # code in lsst throughputs with mods to take out a factor of 3631 to convert
    # the units to closer to "normal" flamba and to write out the un-normalized
    # spectra for use as a SourceSpec.
    sso_B_file = ConfigItem("spectra/B.dat", "B-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_C_file = ConfigItem("spectra/C.dat", "C-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_D_file = ConfigItem("spectra/D.dat", "D-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_Q_file = ConfigItem("spectra/Q.dat", "Q-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_S_file = ConfigItem("spectra/S.dat", "S-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_V_file = ConfigItem("spectra/V.dat", "V-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")
    sso_X_file = ConfigItem("spectra/X.dat", "X-type asteroid from Bus-DeMeo taxonomy (DeMeo et al. 2009)")

    source_mapping = {
                        'sun'    : sun_file,
                        'vega'   : vega_file,
                        'sso::B' : sso_B_file,
                        'sso::C' : sso_C_file,
                        'sso::D' : sso_D_file,
                        'sso::Q' : sso_Q_file,
                        'sso::S' : sso_S_file,
                        'sso::V' : sso_V_file,
                        'sso::X' : sso_X_file
                     }

    sky_brightness_file = ConfigItem("comp/Sky_brightness.dat", "Walker (1987) Sky brightness model")
#    pickles_library_path = ConfigItem("https://ssb.stsci.edu/trds/grid/pickles/dat_uvi/", "Path to Pickles UVILIB spectral library")
    pickles_library_path = ConfigItem('$CDBS_PATH/calspec/pickles', "Path to Pickles UVILIB spectral library")

conf = Conf()

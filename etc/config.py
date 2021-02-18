import os
from astropy.config import ConfigNamespace, ConfigItem

__all__ = ['conf']


class Conf(ConfigNamespace):
    """Configuration parameters."""

    lco_u_file = ConfigItem('comp/lco/SDSS.up.txt', 'LCO SDSS u')
    lco_g_file = ConfigItem('comp/lco/SDSS.gp.txt', 'LCO SDSS g')
    lco_r_file = ConfigItem('comp/lco/SDSS.rp.txt', 'LCO SDSS r')
    lco_i_file = ConfigItem('comp/lco/SDSS.ip.txt', 'LCO SDSS i')
    lco_zs_file = ConfigItem('comp/lco/PSTR-ZS-avg.txt', 'LCO SDSS/PanSTARRS zs')

    lco_c2_file = ConfigItem("comp/lco/LCO_ESA_C2.csv", "LCO ESA C2")
    lco_c3_file = ConfigItem("comp/lco/LCO_ESA_C3.csv", "LCO ESA C3")
    lco_oh_file = ConfigItem("comp/lco/LCO_ESA_OH.csv", "LCO ESA OH")
    lco_cn_file = ConfigItem("comp/lco/LCO_ESA_CN.csv", "LCO ESA CN")
    lco_nh2_file = ConfigItem("comp/lco/LCO_ESA_NH2.csv", "LCO ESA NH2")
    lco_cr_file = ConfigItem("comp/lco/LCO_ESA_CR.csv", "LCO ESA CR")

    lco_U_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-ux.txt', 'LCO Bessell U')
    lco_B_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-bx.txt', 'LCO Bessell B')
    lco_V_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-vx.txt', 'LCO Bessell V')
    lco_R_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-rx.txt', 'LCO Bessell R')
    lco_I_file = ConfigItem('$CDBS_PATH/comp/lco/bssl-ix.txt', 'LCO Bessell I')

    wht_U_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.RGO_U9', 'WHT/PFIP RGO U9')
    wht_B_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_B', 'WHT/PFIP Harris B')
    wht_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_V', 'WHT/PFIP Harris V')
    wht_R_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_R', 'WHT/PFIP Harris R')
    wht_I_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_I', 'WHT/PFIP Harris I')

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

    optics_NaCl_file = ConfigItem("comp/optics/NaCl.dat", "NaCl")
    optics_BAK2_file = ConfigItem("comp/optics/BAK2_glass.dat", "BAK2 glass")
    optics_CaF2_file = ConfigItem("comp/optics/CaF2_glass.dat", "CaF2 glass")
    optics_FK5_file = ConfigItem("comp/optics/FK5_glass.dat", "FK5 crown glass")
    optics_FK58_file = ConfigItem("comp/optics/FK58_glass.dat", "FK58 crown glass")
    optics_LAL7_file = ConfigItem("comp/optics/LAL7_glass.dat", "LAL7 glass")
    optics_UVFS_file = ConfigItem("comp/optics/UVFusedSilica.dat", "UV Fused Silica")

    mapping = { 'u' : lco_u_file,
                'g' : lco_g_file,
                'r' : lco_r_file,
                'i' : lco_i_file,
                'z' : lco_zs_file,
                'zs' : lco_zs_file,
                'C2' : lco_c2_file,
                'C3' : lco_c3_file,
                'OH' : lco_oh_file,
                'CN' : lco_cn_file,
                'NH2': lco_nh2_file,
                'CR' : lco_cr_file,
                'U' : lco_U_file,
                'B' : lco_B_file,
                'V' : lco_V_file,
                'R' : lco_R_file,
                'I' : lco_I_file,
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
    sky_brightness_file = ConfigItem("comp/Sky_brightness.dat", "Walker (1987) Sky brightness model")

conf = Conf()

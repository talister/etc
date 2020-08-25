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

    wht_V_file = ConfigItem('http://svo2.cab.inta-csic.es/theory/fps/getdata.php?format=ascii&id=WHT/PFIP.Har_V', 'WHT/PFIP Harris V')
    # STANDARD STARS
    vega_file = ConfigItem(
        'http://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_010.fits', 'Vega')
    sun_file = ConfigItem(os.path.join('$CDBS_PATH', 'calspec', 'sun_reference_stis_002.fits'), "Solar reference spectrum from https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec.html")

conf = Conf()

from astropy.config import ConfigNamespace, ConfigItem

__all__ = ['conf']


class Conf(ConfigNamespace):
    """Configuration parameters."""

    lco_u_file = ConfigItem('comp/lco/SDSS.up.txt', 'LCO SDSS u')
    lco_g_file = ConfigItem('comp/lco/SDSS.gp.txt', 'LCO SDSS g')
    lco_r_file = ConfigItem('comp/lco/SDSS.rp.txt', 'LCO SDSS r')
    lco_i_file = ConfigItem('comp/lco/SDSS.ip.txt', 'LCO SDSS i')
    lco_zs_file = ConfigItem('comp/lco/PSTR-ZS-avg.txt', 'LCO SDSS/PanSTARRS zs')


conf = Conf()

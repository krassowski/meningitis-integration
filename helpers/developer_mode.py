from jupyter_helpers.notifications import Notifications

Notifications(
    success_audio='/usr/share/sounds/gnome/default/alerts/drip.ogg', time_threshold=3,
    failure_audio='/usr/share/sounds/gnome/default/alerts/bark.ogg',
    integration='GNOME'
)

{% extends 'pbs-base.tpl' %}

{% block environment %}
{% include 'pbs-py3.tpl' %}
{% endblock %}

{% block commands %}
partition-streaming.py -k {{K|default('27')}} -N {{N|default('4')}} -x {{tablesize|default('2e9')}} --stats-interval {{stats_interval|default('10000')}} --stats-dir {{stats_dir}} --pairing-mode {{pairing}} {{samples|join(' ')}}
{% endblock %}

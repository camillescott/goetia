#!/bin/bash -login

#PBS -l walltime={{time|default('02:00:00')}}
#PBS -l nodes={{nodes|default('1')}}:ppn={{ppn|default('1')}}
#PBS -l mem={{mem|default('4gb')}}

#PBS -r n
#PBS -m abe
#PBS -W umask=027
#PBS -N {{name|default('basejob')}}

{% if account is defined %}
#PBS -A {{account}}
{% endif %}
{% if email is defined %}
#PBS -M {{email}}
{% endif %}

{% block environment %}
{% endblock %}

{% block commands %}
{% endblock %}

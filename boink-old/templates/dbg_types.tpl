{# boink/templates/dbg_types.tpl
 # Copyright (C) 2018 Camille Scott
 # All rights reserved.
 #
 # This software may be modified and distributed under the terms
 # of the MIT license.  See the LICENSE file for details.
 #}

{% macro iter_types(storage_types, shifter_types) -%}
{% for Storage_t in storage_types %}
{% for Shifter_t in shifter_types %}
{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}
{{ caller(Storage_t, Shifter_t, tparams, suffix) }}
{% endfor %}
{% endfor %}
{% endmacro %}

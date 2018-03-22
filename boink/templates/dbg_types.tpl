{% macro iter_types(storage_types, shifter_types) -%}
{% for Storage_t in storage_types %}
{% for Shifter_t in shifter_types %}
{% set tparams %}{{Storage_t}},{{Shifter_t}}{% endset %}
{% set suffix %}{{Storage_t}}_{{Shifter_t}}{% endset %}
{{ caller(Storage_t, Shifter_t, tparams, suffix) }}
{% endfor %}
{% endfor %}
{% endmacro %}

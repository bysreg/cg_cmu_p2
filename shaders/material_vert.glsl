varying vec3 norm;
varying vec3 cam_dir;
varying vec3 color;

void main(void)
{
	norm = normalize(gl_NormalMatrix * gl_Normal);
	cam_dir = normalize((gl_ModelViewMatrix*gl_Vertex).xyz);
	//norm = gl_Normal.xyz;
	//norm = normalize((gl_ModelViewProjectionMatrix * vec4(norm, 0.0)).xyz);
	//cam_dir = normalize((gl_ProjectionMatrix * (gl_ModelViewMatrix*gl_Vertex)).xyz);

	gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;
	vec4 col = gl_Color;
	if (all(equal(gl_Color, vec4(1.0)))) col = vec4(0.0);
	gl_FrontColor = col;
	gl_BackColor = col;
}

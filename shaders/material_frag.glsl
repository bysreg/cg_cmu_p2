varying vec3 norm;
varying vec3 cam_dir;
varying vec3 color;

uniform samplerCube u_texture_cube;

void main(void)
{
	vec3 reflect_vector = reflect(cam_dir, norm);
	gl_FragColor = textureCube(u_texture_cube, reflect_vector);
}
